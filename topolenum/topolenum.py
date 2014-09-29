import copy,math,itertools
from sys import stderr
from collections import defaultdict,namedtuple,Counter
from Bio.Phylo.BaseTree import Tree, Clade
from .tree import T

#===============================================================================
# Bugfixes on this branch
#===============================================================================

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
      assemblyobj = copy.deepcopy(assemblyobj)
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
    assemblyobj.recompute()
    assemblyobj.score += self.score
    return assemblyobj


class TreeAssembly(object):
  
  IndexedClade = namedtuple('IndexedClade',['index','clade'])
  
  class KeyPassingDefaultDict(defaultdict):
    def __missing__(self,key):
      self[key] = self.default_factory(key)
      return self[key]
  
  def __init__(self,pwleafdist_histograms,constraint_freq_cutoff,leaves_to_assemble,absolute_freq_cutoff=0.01):
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
                                           key=lambda x: (x.dist,1-x.freq) # sorted on (shortest dist, highest freq)
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
  
  def __deepcopy__(self,memodict={}):
    copy_of_self = type(self).__new__(type(self))
    copy_of_self._nested_set_reprs = [r for r in self._nested_set_reprs]
    copy_of_self._distances_to_root = {k:v for k,v in self._distances_to_root.iteritems()}
    copy_of_self.free_leaves = {fl for fl in self.free_leaves}
    copy_of_self.score = self.score
    copy_of_self._pairs_accounted_for = {p for p in self._pairs_accounted_for}
    copy_of_self.built_clades = [bc for bc in self.built_clades]
    copy_of_self.constraints_idx = [c for c in self.constraints_idx]
    return copy_of_self
  
  def recompute(self,*args):
    if not args or '_distances_to_root' in args:
      self._distances_to_root = {leaf:clade.trace_dist(leaf) for clade in self.built_clades
                                 for leaf in clade.leaf_names}
    if not args or '_pairs_accounted_for' in args:
      self._pairs_accounted_for = {frozenset(pair) for clade in self.built_clades
                                   for pair in itertools.combinations(clade.leaf_names,2)}
    if not args or '_nested_set_reprs' in args:
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
    elif hasattr(extension,'built_clade'):
      new_clade = frozenset({frozenset({extension.built_clade.clade.nested_set_repr(),
                                        extension.new_leaf}),'r'})
      indeces_to_skip = {extension.built_clade.index}
    else:
      new_clade = frozenset({frozenset(c.clade.nested_set_repr()
                                       for c in extension.clades),
                             'r'})
      indeces_to_skip = {c.index for c in extension.clades}
    
    
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
                                 previously_seen,min_score=None):
    joins = self.verify_remaining_proposed_pairs(joins)
    attachments = self.verify_remaining_proposed_pairs(attachments)
    
    for extension_set in (new_pairs,joins,attachments):
      for key,extension in extension_set.items():
        # First filter: has extension been encountered before?
        ext_nested_set_repr = self.as_nested_sets(extension)
        if ext_nested_set_repr in previously_seen:
          extension_set.pop(key)
          continue
        else:
          previously_seen.add(ext_nested_set_repr)
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
  
  def find_extensions(self,previously_seen,min_score=None):
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
    return self.filter_proposed_extensions(new_pairs,joins,attachments,previously_seen,min_score)
    
  def build_extensions(self,new_pairs,joins,attachments):
    # Will need key of pair in new pairs, but not of keys in joins or attachments
    all_ext_to_build = new_pairs.items()+joins.values()+attachments.values()
    updated_assemblies = []
    while all_ext_to_build:
      extension = all_ext_to_build.pop()
      try:
        if not all_ext_to_build:
          updated_assemblies.append(extension.build_extension(self,in_place=True))
        else:
          updated_assemblies.append(extension.build_extension(self))
      except AttributeError:
        if not all_ext_to_build:
          build_in = self
        else:
          build_in = copy.deepcopy(self)
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
        build_in.recompute()
        build_in.score += math.log(pair.freq)
        updated_assemblies.append(build_in)
    return updated_assemblies
  
  def generate_extensions(self,encountered_assemblies,min_score=None):
    new_pairs,joins,attachments = self.find_extensions(encountered_assemblies,min_score)
    if any((new_pairs,joins,attachments)):
      # Return built extensions and any additional information
      return self.build_extensions(new_pairs, joins, attachments)
    else:
      # Somehow inform caller that this assembly is unextendable - it is the caller's responsibility to
      # remove this assembly from the container of active assemblies
      return None

def prepare_assemblies_for_iteration(accepted_assemblies,num_requested_trees,assemblies_workspace,
                                     processing_bite_size,iternum):
  # If we don't yet have a sufficient number of completed assemblies, then we can't filter proposed
  # extensions on min_score, reducing the efficiency of our search. So we need to hurry up and
  # complete at least as many assemblies as have been requested. If we are already working on more
  # that requested number of assemblies, we'll sort them to bring the most likely to get finished
  # first to the front and work this iteration only on as many as were requested to speed up iterations
  # until we have accepted enough trees.
  special_sort1 = lambda x: (len(x.free_leaves),sum((2*c.count_terminals()-1)/2
                            for c in x.built_clades)/x.score)
  special_sort2 = lambda x: sum((2*c.count_terminals()-1)/2 for c in x.built_clades)/x.score
  if len(accepted_assemblies) < num_requested_trees and len(assemblies_workspace) > num_requested_trees:
    # Sorting 1) by number of attached leaves, then 2) by avg score of pairwise distance of all pairs
    # of attached leaves where both are attached under common root.
    return sorted(assemblies_workspace,key=special_sort1)[:num_requested_trees]
  # If sufficient number have been accepted:
  elif len(accepted_assemblies) >= num_requested_trees:
    # On every second iteration sort only on avg score of pw distances disregarding the number of remaining
    # leaves, but take the full set of working assemblies or processing_bite_size many assemblies, whichever
    # one is smaller.
    if iternum % 2 == 0:
      return sorted(assemblies_workspace,key=special_sort1)[:processing_bite_size]
    else:
      # On every sixth iteration sort in avg score order, disregarding the number of remaining leaves, and
      # reverse to let the laggers at the end catch up, but only take half as many assemblies to work on.
      if iternum % 3 == 0:
        return sorted(assemblies_workspace,key=special_sort2,
                      reverse=True)[:processing_bite_size/2]
      # Remaining one third of iterations sort in avg score order and take processing_bite_size many assemblies.
      else:
        return sorted(assemblies_workspace,key=special_sort2)[:processing_bite_size]
  else:
    return [a for a in assemblies_workspace]

def assembly_iteration_new(assemblies_workspace,completed_assemblies,accepted_assemblies,
                           encountered_assemblies,num_requested_trees,processing_bite_size,
                           iternum):
  added_this_iteration = 0
  completed_this_iteration = 0
  discarded_this_iteration = 0
  
  work_on_this_iteration = prepare_assemblies_for_iteration(accepted_assemblies, num_requested_trees,
                                                            assemblies_workspace, processing_bite_size,
                                                            iternum)
  
  print >>stderr,"On iteration",iternum+1," -",len(assemblies_workspace)+len(accepted_assemblies)+len(completed_assemblies),\
                  "TOTAL total assemblies"
  print >>stderr,len(assemblies_workspace)+len(accepted_assemblies),\
                 "total assemblies,",len(accepted_assemblies)+len(completed_assemblies),"completed, working on",\
                 len(work_on_this_iteration),"this iteration, total of",\
                 len(assemblies_workspace),"are incomplete"
  
  assemblies_per_dot = int(math.ceil(float(len(work_on_this_iteration))/50))
  drop_these_from_workspace = []
  for i,assembly in enumerate(work_on_this_iteration):
    print >>stderr,'['+'|'*(i % 10)+' '*(9-(i % 10))+']'+'['+'*'*(i/assemblies_per_dot)+' '*(50-i/assemblies_per_dot)+']'+\
                        '\tcompleted:',completed_this_iteration,'discarded:',discarded_this_iteration,\
                        'new:',added_this_iteration,'\r',
    
    if len(accepted_assemblies) >= num_requested_trees:
      extended_assemblies = assembly.generate_extensions(encountered_assemblies,
                                                         min_score=accepted_assemblies[-1].score)
    else:
      extended_assemblies = assembly.generate_extensions(encountered_assemblies)
    if extended_assemblies is None:
      discarded_this_iteration += 1
      drop_these_from_workspace.append(i)
      continue
    else:
      assert assembly is extended_assemblies[-1]
    
    drop_these = []
    for j,ext_assembly in enumerate(extended_assemblies):
      if ext_assembly.complete:
        completed_this_iteration += 1
        drop_these.append(j)
        if ext_assembly is assembly:
          drop_these_from_workspace.append(i)
        accepted_assemblies.append(ext_assembly)
        accepted_assemblies = sorted(accepted_assemblies,key=lambda x: x.score,reverse=True)
        if len(accepted_assemblies) > num_requested_trees:
          assert len(accepted_assemblies)-1 == num_requested_trees
          completed_assemblies.append(accepted_assemblies.pop())
      elif ext_assembly is not assembly:
        assemblies_workspace.append(ext_assembly)
        added_this_iteration += 1
      else:
        assert j == len(extended_assemblies)-1
#     for j in drop_these[::-1]:
#       extended_assemblies.pop(j)
  
  for i in drop_these_from_workspace[::-1]:
    assemblies_workspace.pop(assemblies_workspace.index(work_on_this_iteration[i]))
  
  print >>stderr,' '*110+'\r',
  if completed_this_iteration:
    print >>stderr,'\t',completed_this_iteration,"were completed this iteration"
  if added_this_iteration:
    print >>stderr,'\t',added_this_iteration,"were added to the workspace this iteration"
  if discarded_this_iteration:
    print >>stderr,'\t',discarded_this_iteration,"were discarded from the workspace this iteration"
  if accepted_assemblies:
    print >>stderr,'\t'+"Best accepted score after this iteration is",accepted_assemblies[0].score
    print >>stderr,'\t'+"Worst accepted score after this iteration is",accepted_assemblies[-1].score
  print >>stderr,'\t',len(encountered_assemblies),"unique assemblies have been encountered"
  
  return accepted_assemblies
    

def assemble_trees(pwleafdist_histograms,leaves_to_assemble,num_requested_trees=1000,
                   constraint_freq_cutoff=0.9,num_iter=10,processing_bite_size=10000,
                   absolute_freq_cutoff=0.01):
  # Init parent of all assemblies
  assemblies_workspace = [TreeAssembly(pwleafdist_histograms,constraint_freq_cutoff,
                                       leaves_to_assemble,absolute_freq_cutoff)]
  accepted_assemblies = []
  completed_assemblies = []
  encountered_assemblies = set()
  iterations = 0
  enough_assemblies_accepted = False
  while assemblies_workspace and iterations < num_iter:
    updated_accepted_assemblies = assembly_iteration_new(assemblies_workspace,completed_assemblies,
                                                         accepted_assemblies,encountered_assemblies,
                                                         num_requested_trees,processing_bite_size,
                                                         iterations)
    if updated_accepted_assemblies is not accepted_assemblies:
      accepted_assemblies = updated_accepted_assemblies
    
    # Did we reach the requested number of assemblies on this iteration?
    if not enough_assemblies_accepted and len(accepted_assemblies) >= num_requested_trees:
      # If yes, we should remove from the workspace assemblies that, unfinished, already have worse
      # scores than the worst accepted assembly
      prev_workspace_length = len(assemblies_workspace)
      assemblies_workspace = filter(lambda x: x.score > accepted_assemblies[-1].score,
                                    assemblies_workspace)
      print >>stderr,"*** Reached",num_requested_trees,"accepted trees, dumping",\
                      prev_workspace_length - len(assemblies_workspace),\
                      "assemblies with scores worse than the worst accepted tree ***"
      enough_assemblies_accepted = True
    iterations += 1
  return assemblies_workspace,accepted_assemblies
