import math,itertools,tempfile,weakref
import cPickle as pickle
from sys import stderr
from collections import defaultdict,namedtuple
from Bio.Phylo.BaseTree import Tree, Clade
from .tree import T as T_BASE


class T(T_BASE):
  
  def __new__(cls,*args,**kwargs):
    pass
  
  def _spawn(self,what_to_spawn,host=None,format=None):
    pass
  
  def _get_spawned(self,item,host=None,format=None):
    pass


LPDF = namedtuple('LeafPairDistanceFrequency',['leaves','dist','freq'])


class ProposedExtension(object):
  
  IndexedPair = namedtuple('IndexedPair',['index','pair'])
  
  def __init__(self,child1,child2):
    if any([isinstance(child1,str),isinstance(child2,str)]):
      # This extension is an attachment of a new leaf to a built clade ...
      # ... so there should be one of each.
      assert isinstance(child1,str) != isinstance(child2,str)
      self.new_leaf = child1 if isinstance(child1,str) else child2 # Figure out ...
      self.built_clade = child2 if isinstance(child1,str) else child1 # .. which is which
      # When a new leaf is attached to a built clade via a new root, the distance between
      # the new leaf and each of the leaves in the built clade will be distance of leaf
      # in existing clade to its current root + 1 to account for the new root
      self.unverified = dict((frozenset({leaf,self.new_leaf}),
                                     self.built_clade.clade.trace_dist(leaf)+1
                              )
                                    for leaf in self.built_clade.clade.leaf_names
                             )
    else: # This extension is the joining of two built clades
      self.clades = sorted([child1,child2],key=lambda x: x.index,reverse=True)
      # When two built clades are joined, the resulting distance between any pair of leaves
      # such that each leaf belongs to a different clade will be
      # sum(distance of each leaf to its current root) + 1 to account for the new root.
      self.unverified = dict((frozenset({leafpair[0][0],leafpair[1][0]}),
                                              leafpair[0][1]+leafpair[1][1]+1
                                     )
                                    for leafpair
                                    in itertools.product(*([(leaf,
                                                             clade.clade.trace_dist(leaf)
                                                             )
                                                            for leaf in clade.clade.leaf_names
                                                            ]
                                                           for clade in self.clades
                                                           )
                                                         )
                                    )
    self.consistent = {}
    self.inconsistent = {}
    self.verified = set()
    self.score = 0.0
  
  def check_pair(self,pair,i):
    if pair.leaves in self.consistent:
      # If this leaf pair has already been added to consistent, then ...
      # ... it should not be in unverified any more ...
      assert pair.leaves not in self.unverified
      # ... and this new distance should be different from the consistent one ...
      assert pair.dist != self.consistent[pair.leaves].pair.dist
      self.inconsistent[i] = pair # ... so it goes into inconsistent
    else: # If it hasn't been added to consistent ...
      if pair.dist == self.unverified[pair.leaves]: # ... and its distance matches expected
        self.consistent[pair.leaves] = self.IndexedPair(i,pair) # it goes into consistent
        # and is "verified", so pop it from unverified and add to verified
        self.unverified.pop(pair.leaves)
        self.verified.add(pair.leaves)
        self.score += math.log(pair.freq)
      else: # ... and its distance doesn't match expected
        self.inconsistent[i] = pair # it goes into inconsistent ...
        # ... and it remains "unverified", so don't pop it from unverified
  
  def build_extension(self,assemblyobj,in_place=False):
    if not in_place:
      assemblyobj = assemblyobj.copy()
    assert not self.unverified
    
    # Assemble for removal from constraints_idx indeces of consistent and inconsistent pairdists
    drop_these = [pair.index for pair in self.consistent.values()]+self.inconsistent.keys()
    
    # Make sure to use clade(s) from this assemblyobj for construction, not the ones from
    # the original that are stored in self.clades
    if hasattr(self,'new_leaf'):
      built_clade = assemblyobj.built_clades.pop(self.built_clade.index)
      assemblyobj.free_leaves.remove(self.new_leaf)
      assert set(self.built_clade.clade.leaf_names) == set(built_clade.leaf_names)
      new_clades_attr = [built_clade.wrapped,Clade(name=self.new_leaf)]
      # One more thing to do if this extension is an attachment of a new leaf:
      # Add for removal constraint pairdists where the new leaf has distance 1 with any
      # other leaf, regardless of whether the other leaf is involved in this extension.
      # All such constraint distances are inconsistent with this extension. This additional
      # check is necessary because of the special way new pair extensions are handled - with
      # no questions asked.
      drop_these.extend([i for i in assemblyobj.constraints_idx
                         if self.new_leaf in assemblyobj.constraints_master[i].leaves and
                                             assemblyobj.constraints_master[i].dist == 1])
    else:
      assert all(set(c.clade.leaf_names)==set(assemblyobj.built_clades[c.index].leaf_names)
                                                            for c in self.clades)
      new_clades_attr = [assemblyobj.built_clades.pop(c.index).wrapped for c in self.clades]
    
    assemblyobj.constraints_idx = filter(lambda x: x not in drop_these,
                                         assemblyobj.constraints_idx)
    
    assemblyobj.built_clades.append(T(Clade(clades=new_clades_attr)))
    assemblyobj.recompute(extension=self)
    assemblyobj.score += self.score
    return assemblyobj


class TreeAssembly(object):
  
  IndexedClade = namedtuple('IndexedClade',['index','clade'])
  
  class KeyPassingDefaultDict(defaultdict):
    def __missing__(self,key):
      self[key] = self.default_factory(key)
      return self[key]
  
  def __init__(self,pwleafdist_histograms,constraint_freq_cutoff,leaves_to_assemble,
               absolute_freq_cutoff=0.01):
    #===========================================================================
    # The data attributes below will are set on the class, not in instances,
    # meaning they will be shared between all instances of this class, saving
    # lots of space as the number of assemblies grows.
    # *** The are INVARIANTS and should not be changed by instances!!! ***
    #===========================================================================
    
    # Make a sorted master reference tuple of pw distance freq constraints, ...
    type(self).constraints_master = tuple(sorted([
                                            # by creating LPDF tuples ...
                                            LPDF(pair_histogram[0],score[0],score[1])
                                            # for every histogram of pw distances for a leaf pair ...
                                            for pair_histogram in
                                            # present in a built-in-place subset of
                                            # the full histogram for that pair ...
        [(leafpair,[dist for i,dist in enumerate(distances) if
                    # of distances ordered by freq such that the sum of their freqs
                    # is less than the requested frequency cutoff, ...
                    sum(d[1] for d in distances[:i]) < constraint_freq_cutoff]
          ) for leafpair,distances in pwleafdist_histograms]
                                            for score in pair_histogram[1]],
                                           # sorted on (shortest dist, highest freq)
                                           key=lambda x: (x.dist,1-x.freq)
                                                 )
                                          )
    
    # Unlike constraints_master, which guarantees invariance by being recursively immutable,
    # this is a mutable dict. Is there a way to force immutability?
    type(self).pwdist_histograms_dict = {leafpair:dict(dist_histogram)
                                         for leafpair,dist_histogram in pwleafdist_histograms}
    
    # Again, a mutable variable, though I would like it to be immutable is possible
    type(self).abs_cutoff = absolute_freq_cutoff
    
    #===========================================================================
    # END of class attributes
    #===========================================================================
    
    self.built_clades = []
    self.free_leaves = set(leaves_to_assemble)
    self.constraints_idx = range(len(self.constraints_master))
    self.score = 0.0
  
  def __getstate__(self):
    pass
  
  def __setstate__(self,state):
    pass
  def copy(self):
    copy_of_self = type(self).__new__(type(self))
    copy_of_self._nested_set_reprs = [r for r in self._nested_set_reprs]
    copy_of_self._distances_to_root = {k:v for k,v in self._distances_to_root.iteritems()}
    copy_of_self.free_leaves = {fl for fl in self.free_leaves}
    copy_of_self.score = self.score
    copy_of_self._pairs_accounted_for = {p for p in self._pairs_accounted_for}
    copy_of_self.built_clades = [bc for bc in self.built_clades]
    copy_of_self.constraints_idx = [c for c in self.constraints_idx]
    return copy_of_self
  
  def recompute(self,*args,**kwargs):
    if 'extension' in kwargs:
      extension = kwargs['extension']
      if type(extension).__name__ == 'LeafPairDistanceFrequency':
        for leaf in extension.leaves:
          self._distances_to_root[leaf] = 1
        self._pairs_accounted_for.add(extension.leaves)
        self._nested_set_reprs.append(frozenset({frozenset(extension.leaves),'r'}))
      else:
        self._distances_to_root = extension.distances_to_root
        self._pairs_accounted_for = extension.pairs_accounted_for
        self._nested_set_reprs = extension.nested_set_reprs
    else:
      if '_distances_to_root' in args:
        self._distances_to_root = {leaf:clade.trace_dist(leaf) for clade in self.built_clades
                                   for leaf in clade.leaf_names}
      if '_pairs_accounted_for' in args:
        self._pairs_accounted_for = {frozenset(pair) for clade in self.built_clades
                                     for pair in itertools.combinations(clade.leaf_names,2)}
      if '_nested_set_reprs' in args:
        self._nested_set_reprs = [frozenset({c.nested_set_repr(),'r'}) for c in self.built_clades]
  
  def _property_getter(self,property):
    try:
      return getattr(self,property)
    except AttributeError:
      self.recompute(property)
      return getattr(self,property)
  
  @property
  def current_clades_as_nested_sets(self):
    return self._property_getter('_nested_set_reprs')
  
  @property
  def distances_to_root(self):
    return self._property_getter('_distances_to_root')
  
  @property
  def pairs_accounted_for(self):
    return self._property_getter('_pairs_accounted_for')
  
  @property
  def complete(self):
    return len(self.built_clades) == 1 and not self.free_leaves
  
  def verify_remaining_proposed_pairs(self,extensions):
    for key,ext in extensions.items():
      for pair,dist in ext.unverified.items():
        # Looking up pair histogram separately in case lookup throws KeyError,
        # since we don't want to catch that one
        pair_histogram = self.pwdist_histograms_dict[pair]
        try:
          pair_freq = pair_histogram[dist]
        except KeyError:
          # If the resulting pairdist is not in the histogram, it means it wasn't observed
          # at all, so it's frequency is 0.
          pair_freq = 0.0
        if pair_freq < self.abs_cutoff:
          extensions.pop(key)
          break
        else:
          ext.score += math.log(pair_freq)
          ext.unverified.pop(pair)
          ext.verified.add(pair)
      else:
        # If all unverified pairs check out, make sure this extension is not passing
        # this filter entirely on the strength of pairs we just verified
        if not ext.consistent:
          extensions.pop(key)
    return extensions
  
  def as_nested_sets(self,extension):
    # Clades are represented by {clade.nested_set_repr(),'r'} (r for root) to indicate
    # their free-standing nature. This way {{clade1,'r'},{clade2,'r'}} represents two
    # free-standing clades, distinct from {clade1,clade2}, representing a single
    # free-standing clade (or tree) with two sub-clades.
    if hasattr(extension,'freq'):
      new_clade = frozenset({frozenset(extension.leaves),'r'})
      indeces_to_skip = []
    else:
      if hasattr(extension,'built_clade'):
        new_clade = frozenset({frozenset({extension.built_clade.clade.nested_set_repr(),
                                          extension.new_leaf}),'r'})
        indeces_to_skip = {extension.built_clade.index}
      else:
        new_clade = frozenset({frozenset(c.clade.nested_set_repr()
                                         for c in extension.clades),
                               'r'})
        indeces_to_skip = {c.index for c in extension.clades}
      extension.nested_set_reprs = [c for i,c in enumerate(self.current_clades_as_nested_sets)
                                    if i not in indeces_to_skip]
      extension.nested_set_reprs.append(new_clade)
    
    old_clades = frozenset(c for i,c in enumerate(self.current_clades_as_nested_sets)
                           if i not in indeces_to_skip)
    return old_clades|frozenset({new_clade})
  
  def best_case_with_extension(self,extension):
    try:
      updated_distances_to_root = {leaf:(self.distances_to_root[leaf]+1 if leaf in
                                         self.distances_to_root else 1)
                               for pair in itertools.chain(extension.consistent,
                                                           extension.verified)
                               for leaf in pair}
      for leaf in self.distances_to_root:
        if leaf not in updated_distances_to_root:
          updated_distances_to_root[leaf] = self.distances_to_root[leaf]
      
      updated_pairs_accounted_for = {pair for pair in
                                     itertools.chain(self.pairs_accounted_for,
                                                     extension.verified)}
      extension.distances_to_root = updated_distances_to_root
      extension.pairs_accounted_for = updated_pairs_accounted_for
    except AttributeError:
      updated_distances_to_root = dict(self.distances_to_root.iteritems())
      for leaf in extension.leaves:
        updated_distances_to_root[leaf] = 1
      updated_pairs_accounted_for = {pair for pair in self.pairs_accounted_for}
      updated_pairs_accounted_for.add(extension.leaves)
    
    try:
      best_possible_final_score = self.score + extension.score
    except AttributeError:
      best_possible_final_score = self.score + math.log(extension.freq)
    for pair,histogram in self.pwdist_histograms_dict.items():
      if pair not in updated_pairs_accounted_for:
        min_dist = sum(updated_distances_to_root[leaf] if leaf in updated_distances_to_root else
                       0 for leaf in pair)+1
        acceptable_dist_freqs = [freq for dist,freq in histogram.items() if dist >= min_dist]
        if acceptable_dist_freqs:
          best_possible_final_score += math.log(max(acceptable_dist_freqs))
        else:
          return None
    return best_possible_final_score
  
  def filter_proposed_extensions(self,new_pairs,joins,attachments,
                                 encountered,min_score=None):
    joins = self.verify_remaining_proposed_pairs(joins)
    attachments = self.verify_remaining_proposed_pairs(attachments)
    
    for extension_set in (new_pairs,joins,attachments):
      for key,extension in extension_set.items():
        # First filter: has extension been encountered before?
        ext_nested_set_repr = self.as_nested_sets(extension)
        if ext_nested_set_repr in encountered:
          extension_set.pop(key)
          continue
        else:
          encountered.add(ext_nested_set_repr)
        # Second filter: is score with extension already worse than min_score?
        if min_score is not None:
          try: # Will fail is item is a new pair - forgiveness faster than permission
            if extension.score + self.score < min_score:
              extension_set.pop(key)
              continue
          except AttributeError: # Must be a new_pair:
            if math.log(extension.freq) + self.score < min_score:
              extension_set.pop(key)
              continue
        # Third filter: is there a way to extend the extension all the way to a full assembly?
        # Is the upper limit on best score for that assembly already worse than min_score?
        best_case = self.best_case_with_extension(extension)
        if best_case is None or best_case < min_score:
          extension_set.pop(key)
          continue
    
    return new_pairs,joins,attachments
  
  def find_extensions(self,encountered,min_score=None):
    new_pairs = {}
    joins = self.KeyPassingDefaultDict(lambda key: ProposedExtension(*key))
    attachments = self.KeyPassingDefaultDict(lambda key: ProposedExtension(*key))
    already_connected = {frozenset(clade.leaf_names):i for i,clade in
                                                       enumerate(self.built_clades)}
    # All pairwise intersections should be empty
    assert not any(frozenset.intersection(*leafsetpair)
                   for leafsetpair in itertools.combinations(already_connected.iterkeys(),2))
    already_connected_splat = frozenset({leaf for leafset in already_connected
                                              for leaf in leafset})
    ac_leafdict = {leaf:self.IndexedClade(i,self.built_clades[i])
                   for leafset,i in already_connected.iteritems() for leaf in leafset}
    for i in self.constraints_idx:
      pair = self.constraints_master[i]
      if pair.dist == 1:
        # Pairs with distance 1 are added w/o questions. If continue with this path,
        # later we will make sure to remove from consideration all pairs that conflict this.
        new_pairs[i] = pair
      elif not pair.leaves & already_connected_splat:
        # If a pair has distance > 1 and neither leaf in pair has already been added to a
        # clade, then we can't do anything with it, so we silently skip it
        continue
      else:
        if pair.leaves <= already_connected_splat: # If both leaves are already in built clades
          if any(pair.leaves <= built_leafset for built_leafset in already_connected.iterkeys()):
            # If one built clade already contains both leaves, there is nothing to do either
            continue
          else:
            # If separate built clades contain one leaf each, then the pair goes into
            # the corresponding join's ProposedExtension object
            joins[frozenset(ac_leafdict[leaf] for leaf in pair.leaves)].check_pair(pair,i)
        else:
          # (pair.leaves & already_connected_splat) and (not pair.leaves <= already_connected_splat)
          # implies one leaf in pair is already contained in a built clade, the other isn't
          for leaf in pair.leaves: # Identify the new leaf and the built clade
            try:
              clade_of_attached_leaf = ac_leafdict[leaf]
              attached_leaf = leaf
            except KeyError as e:
              new_leaf = leaf
          # And put the pair into the corresponding attachment's ProposedExtension object
          attachments[frozenset({clade_of_attached_leaf,new_leaf})].check_pair(pair,i)
    return self.filter_proposed_extensions(new_pairs,joins,attachments,encountered,min_score)
    
  def build_extensions(self,new_pairs,joins,attachments):
    # Will need key of pair in new pairs, but not of keys in joins or attachments
    all_ext_to_build = new_pairs.items()+joins.values()+attachments.values()
    extended_assemblies = []
    while all_ext_to_build:
      extension = all_ext_to_build.pop()
      try:
        if not all_ext_to_build:
          extended_assemblies.append(extension.build_extension(self,in_place=True))
        else:
          extended_assemblies.append(extension.build_extension(self))
      except AttributeError:
        if not all_ext_to_build:
          build_in = self
        else:
          build_in = self.copy()
        idx_of_pair,pair = extension # Now we can get the key (index of pair in constraints_idx)
        
        # Select for dropping all pairs with distance 1 and one member of pair - they can't
        # have distance 1 with anyone except each other
        drop_these = [i for i in self.constraints_idx if self.constraints_master[i].dist == 1 and
                                              self.constraints_master[i].leaves & pair.leaves and
                                            not self.constraints_master[i].leaves == pair.leaves]
        # Select for dropping all pairs of these two leaves with distance > 1
        drop_these.extend(i for i in self.constraints_idx if self.constraints_master[i].dist > 1
                                            and self.constraints_master[i].leaves == pair.leaves)
        drop_these.append(idx_of_pair) # Finally, select for dropping this pair
        # Drop selected pairs from constraints_idx
        build_in.constraints_idx = filter(lambda x: x not in drop_these,build_in.constraints_idx)
        
        # Remove leaves in this pair from free_leaves
        for leaf in pair.leaves:
          build_in.free_leaves.remove(leaf) # Pop each leaf in pair from free_leaves
        
        # Build new clade and update the score
        build_in.built_clades.append(T(Clade(clades=[Clade(name=leaf) for leaf in pair.leaves])))
        build_in.recompute(extension=pair)
        build_in.score += math.log(pair.freq)
        extended_assemblies.append(build_in)
    return extended_assemblies
  
  def generate_extensions(self,encountered_assemblies,min_score=None):
    new_pairs,joins,attachments = self.find_extensions(encountered_assemblies,min_score)
    if any((new_pairs,joins,attachments)):
      # Return built extensions and any additional information
      return self.build_extensions(new_pairs, joins, attachments)
    else:
      # Somehow inform caller that this assembly is unextendable - it is the caller's responsibility to
      # remove this assembly from the container of active assemblies
      return None


class FIFOfile(object):
  def __init__(self,name='use_tempfile',mode='b',wbuffering=0,rbuffering=0,delete=True):
    if name == 'use_tempfile':
      self._wh = tempfile.NamedTemporaryFile('w'+mode,bufsize=wbuffering,dir='.',delete=delete)
      self.name = self._wh.name
    else:
      self.name = name
      self._wh = open(self.name,'w'+mode,wbuffering)
    self._rh = open(self.name,'r'+mode,rbuffering)
  
  def read(self):
    try:
      result = pickle.load(self._rh)
    except EOFError:
      try:
        self._rh.readline()
        result = pickle.load(self._rh)
      except EOFError:
        return None
    return result
  
  def write(self,item):
    pickle.dump(item,self._wh,pickle.HIGHEST_PROTOCOL)
  
  def close(self):
    self._rh.close()
    self._wh.close()
  


class AssemblyWorkspace(object):
  def __init__(self,pwleafdist_histograms,constraint_freq_cutoff,leaves_to_assemble,
               absolute_freq_cutoff=0.01,num_requested_trees=1000,max_workspace_size=10000,
               tmpfilename=None):
    self.workspace = [TreeAssembly(pwleafdist_histograms,constraint_freq_cutoff,
                                   leaves_to_assemble,absolute_freq_cutoff)]
    self.accepted_assemblies = []
    self.rejected_assemblies = []
    self.encountered_assemblies = set()
    
    self.num_requested_trees = num_requested_trees
    self.reached_num_requested_trees = False
    self.curr_min_score = None
    
    self.iternum = 0
    
    self.max_workspace_size = max_workspace_size
    
    self.filename = tmpfilename
  
  @property
  def _overflow(self):
    try:
      return self._overflowFIFO
    except AttributeError:
      if self.filename is None:
        self._overflowFIFO = FIFOfile('use_tempfile')
      else:
        self._overflowFIFO = FIFOfile(self.filename)
      return self._overflowFIFO
  
  def update_workspace(self,new_assemblies=None):
    while len(self.workspace) < self.max_workspace_size:
      popped = self._overflow.read()
      if popped is None:
        break
      else:
        if popped.score > self.curr_min_score:
          self.workspace.append(popped)
    
    if new_assemblies is not None:
      if len(self.workspace) < self.max_workspace_size:
        while len(self.workspace) < self.max_workspace_size:
          try:
            self.workspace.append(new_assemblies.pop(0))
          except IndexError:
            break
      
      if new_assemblies:
        while new_assemblies:
          self._overflow.write(new_assemblies.pop(0))
  
  def process_extended_assembly(self,assembly):
    if assembly.complete:
      if len(self.accepted_assemblies) < self.num_requested_trees\
                                      or assembly.score > self.curr_min_score:
        self.accepted_assemblies.append(assembly)
        self.accepted_assemblies = sorted(self.accepted_assemblies,
                                          key=lambda x: x.score,reverse=True)
        while len(self.accepted_assemblies) > self.num_requested_trees:
          self.rejected_assemblies.append(self.accepted_assemblies.pop())
        self.curr_min_score = self.accepted_assemblies[-1].score
      else:
        self.rejected_assemblies.append(assembly)
      return
    else:
      return assembly
  
  def iterate(self):
    drop_from_workspace_idx = []
    workspace_this_iter = [assembly for assembly in self.workspace]
    for i,assembly in enumerate(workspace_this_iter):
      extended_assemblies = assembly.generate_extensions(self.encountered_assemblies,
                                                         self.curr_min_score)
      if extended_assemblies is None:
        drop_from_workspace_idx.append(i)
        continue
      else:
        assert assembly is extended_assemblies[-1]
        if self.process_extended_assembly(extended_assemblies.pop()) is None:
          drop_from_workspace_idx.append(i)
      self.update_workspace([asbly for asbly in extended_assemblies
                             if self.process_extended_assembly(asbly) is not None])
    for i in drop_from_workspace_idx[::-1]:
      self.workspace.pop(i)
    if not self.reached_num_requested_trees:
      if len(self.accepted_assemblies) == self.num_requested_trees:
        self.workspace = [asbly for asbly in self.workspace if asbly.score > self.curr_min_score]
    self.update_workspace()
    self.iternum += 1
