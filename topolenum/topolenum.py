import copy,math,itertools
from sys import stderr
from collections import defaultdict,namedtuple
from Bio.Phylo.BaseTree import Tree, Clade
from .tree import T

def check_against_histograms(extensions,pwhistdict,freq_cutoff=0.01):
    good_extensions = {}
    if all(any(isinstance(f,str) for f in ext) for ext in extensions):
        for att,pairs in extensions.items():
            already_passed = [p[0] for p in pairs]
            distance_pieces = [{ln:len(c.get_path(ln)) for ln in T(c).leaf_names} for c in att if not isinstance(c,str)].pop()
            new_leaf = [c for c in att if isinstance(c,str)].pop()
            distances = {frozenset({f,new_leaf}):d+1 for f,d in distance_pieces.items()}
            if all((d in pwhistdict[p] and pwhistdict[p][d] > freq_cutoff) for p,d in distances.items()):
                good_extensions[att] = pairs
    else:
        for j,pairs in extensions.items():
            already_passed = [p[0] for p in pairs]
            distance_pieces = [{ln:len(c.get_path(ln)) for ln in T(c).leaf_names} for c in j]
            distances = {frozenset({f1,f2}):distance_pieces[0][f1]+distance_pieces[1][f2]+1
                         for f1 in distance_pieces[0] for f2 in distance_pieces[1] if frozenset({f1,f2}) not in already_passed}
            if all((d in pwhistdict[p] and pwhistdict[p][d] > freq_cutoff) for p,d in distances.items()):
                good_extensions[j] = pairs
    return good_extensions

def find_legal_extensions(built_clades,available_pairs):
    new_pairs = []
    joins = defaultdict(list)
    attachments = defaultdict(list)
    already_connected = {frozenset({ln for ln in T(c).leaf_names}):c for c in built_clades}
    already_connected_flat = frozenset({ln for k in already_connected for ln in k})
    for pair in available_pairs:
      if pair[1][0] == 1:
        new_pairs.append(pair)
      elif pair[0] & already_connected_flat:
          if pair[0] < already_connected_flat:
              if any(pair[0] <= acs for acs in already_connected.keys()):
                  continue
              else:
                  curr_hosts = {f:[c for ls,c in already_connected.items() if f in ls].pop() for f in pair[0]}
                  assert len(curr_hosts) == 2
                  if pair[1][0] == sum(len(c.get_path(f)) for f,c in curr_hosts.items())+1:
                      joins[frozenset({c for c in curr_hosts.values()})].append(pair)
          else:
              curr_host = {f:[c for ls,c in already_connected.items() if f in ls].pop() for f in pair[0] if f in already_connected_flat}
              assert len(curr_host) == 1,str(curr_host)+' ** '+str(pair)
              new_leaf = [f for f in pair[0] if f not in curr_host.keys()].pop()
              if pair[1][0] == sum(len(c.get_path(f)) for f,c in curr_host.items())+1:
                  attachments[frozenset({c for c in curr_host.values()}.union([new_leaf]))].append(pair)
    return new_pairs,joins,attachments
          
def remove_inconsistent_pairs(extension,available_pair_idxs,pairs_master,extension_type=None):
    remove_these = []
    if extension_type == 'new_pair':
        remove_these.extend(filter(lambda x: pairs_master[x][1][0] == 1 and pairs_master[x][0] & extension and not pairs_master[x][0] == extension,available_pair_idxs))
        remove_these.extend(filter(lambda x: pairs_master[x][1][0] > 1 and pairs_master[x][0] == extension,available_pair_idxs))
    else:
        leafsets = [[c] if isinstance(c,str) else T(c).leaf_names for c in extension]
        compare_against = {frozenset({f1,f2}) for f1 in leafsets[0] for f2 in leafsets[1]}
        tmp = [c for c in extension if isinstance(c,str)]
        new_leaf = tmp.pop() if tmp else False
        for i in available_pair_idxs:
            if pairs_master[i][0] in compare_against:
                remove_these.append(i)
            elif new_leaf:
                if new_leaf in pairs_master[i][0] and pairs_master[i][1][0] == 1:
                    remove_these.append(i)
    
    return filter(lambda x: x not in remove_these,available_pair_idxs)

def look_ahead(clades,pwhistdict,curr_partial_score):
    scored_pairs = []
    partial_distances = {}
    for c in clades:
        tc = T(c)
        leafnames = tc.leaf_names
        for i,ln1 in enumerate(leafnames):
            partial_distances[ln1] = len(tc.get_path(ln1))
            for ln2 in leafnames[i+1:]:
                scored_pairs.append(frozenset({ln1,ln2}))
    best_possible_full_score = curr_partial_score
    for k,v in pwhistdict.items():
        if k not in scored_pairs:
            min_dist = sum(partial_distances[f] if f in partial_distances else 0 for f in k)+1
            if [freq for dist,freq in v.items() if dist >= min_dist]:
                best_possible_full_score += math.log(max(freq for dist,freq in v.items() if dist >= min_dist))
            else:
                return None
    return best_possible_full_score

def build_extension(extension,assembly,pairs_master,pwhistdict,extension_type=None):
    if extension_type == 'new_pair':
        assembly[2].pop(assembly[2].index(pairs_master.index(extension)))
        for f in extension[0]:
            assembly[1].pop(assembly[1].index(f))
        assembly[4] += math.log(extension[1][1])
        assembly[0].append(Clade(clades=[Clade(name=f) for f in extension[0]]))
    else:
        for p in extension[1]:
            assembly[2].pop(assembly[2].index(pairs_master.index(p)))
        identified_pairs = dict(extension[1])
        new_leaf = False
        existing_clade_copies = []
        for piece in extension[0]:
            if isinstance(piece,str):
                new_leaf = piece
                assembly[1].pop(assembly[1].index(piece))
            else:
                candidates_for_copy_of_existing_clade = [c for c in assembly[0] if set(T(c).leaf_names) == set(T(piece).leaf_names)]
                assert len(candidates_for_copy_of_existing_clade) == 1
                existing_clade_copies.append(candidates_for_copy_of_existing_clade.pop())
        if new_leaf:
            assert len(existing_clade_copies) == 1
            existing_clade = existing_clade_copies.pop()
            for f in T(existing_clade).leaf_names:
                if frozenset({new_leaf,f}) in identified_pairs:
                    assembly[4] += math.log(identified_pairs[frozenset({new_leaf,f})][1])
                else:
                    assembly[4] += math.log(pwhistdict[frozenset({new_leaf,f})][len(existing_clade.get_path(f))+1])
            assembly[0].pop(assembly[0].index(existing_clade))
            assembly[0].append(Clade(clades=[existing_clade,Clade(name=new_leaf)]))
        else:
            assert len(existing_clade_copies) == 2
            for f1 in T(existing_clade_copies[0]).leaf_names:
                for f2 in T(existing_clade_copies[1]).leaf_names:
                    if frozenset({f1,f2}) in identified_pairs:
                        assembly[4] += math.log(identified_pairs[frozenset({f1,f2})][1])
                    else:
                        assembly[4] += math.log(pwhistdict[frozenset({f1,f2})][len(existing_clade_copies[0].get_path(f1))+len(existing_clade_copies[1].get_path(f2))+1])
            for c in existing_clade_copies:
                assembly[0].pop(assembly[0].index(c))
            assembly[0].append(Clade(clades=existing_clade_copies))
#     assembly[6] = {frozenset({frozenset(n.leaf_names) for n in T(c).get_nonterminals()}) for c in assembly[0]}
    assembly[2] = remove_inconsistent_pairs(extension[0],assembly[2],pairs_master,extension_type)

def build_all_extensions(new_pairs,joins,attachments,assembly,pairs_master,pwhistdict):
    build_here = [assembly]
    build_here.extend(copy.deepcopy(assembly) for i in xrange(sum(len(coll) for coll in [new_pairs,joins,attachments])-1))
    bld_idx = 0
    for np in new_pairs:
        build_extension(np,build_here[bld_idx],pairs_master,pwhistdict,'new_pair')
        bld_idx += 1
    for ext,pairs in joins.items()+attachments.items():
        build_extension((ext,pairs),build_here[bld_idx],pairs_master,pwhistdict)
        bld_idx += 1
    build_here.pop(0)
    return build_here

def assembly_iteration(assemblies,pairs_master,pwhistdict,accepted_assemblies,encountered_leaf_subsets,num_requested_trees,processing_bite_size,iternum):
    new_this_iteration = []
    accepted_this_iteration = []
    discarded_this_iteration = {'u':[],'e':[],'cs':[],'ps':[]}
    assemblies_this_iter = [(i,asb) for i,asb in enumerate(assemblies) if (asb[1] or len(asb[0]) > 1)]
    if len(accepted_assemblies) < num_requested_trees and len(assemblies_this_iter) > num_requested_trees:
#         assemblies_this_iter = sorted(assemblies_this_iter,key=lambda x: (len(x[1][1]),x[1][4]))[:num_requested_trees]
        assemblies_this_iter = sorted(assemblies_this_iter,key=lambda x: (len(x[1][1]),sum((c.count_terminals()+(c.count_terminals()-1))/2
                                                                             for c in x[1][0])/x[1][4]))[:num_requested_trees]
    elif len(accepted_assemblies) >= num_requested_trees:
#         assemblies_this_iter = sorted(assemblies_this_iter,key=lambda x: (len(x[1][1]),x[1][4]))[:processing_bite_size]
        if iternum % 2 == 0:
            assemblies_this_iter = sorted(assemblies_this_iter,key=lambda x: (len(x[1][1]),sum((c.count_terminals()+(c.count_terminals()-1))/2
                                                                                               for c in x[1][0])/x[1][4]))[:processing_bite_size]
        else:
            if iternum % 3 == 0:
                assemblies_this_iter = sorted(assemblies_this_iter,key=lambda x: sum((c.count_terminals()+(c.count_terminals()-1))/2
                                                                                     for c in x[1][0])/x[1][4],reverse=True)[:processing_bite_size/2]
            else:
                assemblies_this_iter = sorted(assemblies_this_iter,key=lambda x: sum((c.count_terminals()+(c.count_terminals()-1))/2
                                                                                       for c in x[1][0])/x[1][4])[:processing_bite_size]
    print >>stderr,len(assemblies),"tracked trees,",len(accepted_assemblies),"accepted, working on",\
                    len(assemblies_this_iter),"this iteration, total of",len(assemblies)-len(accepted_assemblies),"are incomplete"
    assemblies_per_dot = int(math.ceil(float(len(assemblies_this_iter))/50))
    for i,assembly in assemblies_this_iter:
        print >>stderr,'['+'|'*(assemblies_this_iter.index((i,assembly)) % 10)+' '*(9-(assemblies_this_iter.index((i,assembly)) % 10))+']'+\
                        '['+'*'*(assemblies_this_iter.index((i,assembly))/assemblies_per_dot)+' '*(50-assemblies_this_iter.index((i,assembly))/assemblies_per_dot)+']'+\
                        '\taccepted:',len(accepted_this_iteration),'discarded:',sum(len(v) for v in discarded_this_iteration.values()),\
                        "new:",len(new_this_iteration),'\r',
        new_pairs,joins,attachments = find_legal_extensions(assembly[0],[pairs_master[i] for i in assembly[2]])
        good_joins = check_against_histograms(joins,pwhistdict)
        good_attachments = check_against_histograms(attachments,pwhistdict)
        if not new_pairs and not good_joins and not good_attachments:
            discarded_this_iteration['e'].append(assemblies.pop(assemblies.index(assembly)))
            continue
        new_assemblies = build_all_extensions(new_pairs,good_joins,good_attachments,assembly,pairs_master,pwhistdict)
        still_in_play = True
        if len(accepted_assemblies) == num_requested_trees and assembly[4] < accepted_assemblies[-1][4]:
            discarded_this_iteration['cs'].append(assemblies.pop(assemblies.index(assembly)))
            encountered_leaf_subsets.add(frozenset({frozenset({frozenset(n.leaf_names) for n in T(c).get_nonterminals()}) for c in assembly[0]}))
            still_in_play = False
        else:
            la = look_ahead(assembly[0],pwhistdict,assembly[4])
            if la is None:
                discarded_this_iteration['e'].append(assemblies.pop(assemblies.index(assembly)))
                encountered_leaf_subsets.add(frozenset({frozenset({frozenset(n.leaf_names) for n in T(c).get_nonterminals()}) for c in assembly[0]}))
                still_in_play
            elif len(accepted_assemblies) == num_requested_trees and la < accepted_assemblies[-1][4]:
                discarded_this_iteration['ps'].append(assemblies.pop(assemblies.index(assembly)))
                encountered_leaf_subsets.add(frozenset({frozenset({frozenset(n.leaf_names) for n in T(c).get_nonterminals()}) for c in assembly[0]}))
                still_in_play = False
            else:
                leafsets = frozenset({frozenset({frozenset(n.leaf_names) for n in T(c).get_nonterminals()}) for c in assembly[0]})
                if leafsets in encountered_leaf_subsets:
                    discarded_this_iteration['u'].append(assemblies.pop(assemblies.index(assembly)))
                    still_in_play = False
                else:
                    encountered_leaf_subsets.add(leafsets)
    #             for other_assembly in assemblies:
#                 if other_assembly is not assembly:
#                     if assembly[6] == other_assembly[6]:
#                         assert set(assembly[2]) == set(other_assembly[2])
#                         discarded_this_iteration['u'].append(assemblies.pop(assemblies.index(assembly)))
#                         still_in_play = False
#                         break
#             else:
#                 if assembly[6] in encountered_leaf_subsets:
#                     discarded_this_iteration['u'].append(assemblies.pop(assemblies.index(assembly)))
#                     still_in_play = False
        rejected_new_assemblies = []
        for new_assembly in new_assemblies:
            if len(accepted_assemblies) == num_requested_trees and new_assembly[4] < accepted_assemblies[-1][4]:
                rejected_new_assemblies.append(new_assembly)
                encountered_leaf_subsets.add(frozenset({frozenset({frozenset(n.leaf_names) for n in T(c).get_nonterminals()}) for c in new_assembly[0]}))
            else:
                la = look_ahead(new_assembly[0],pwhistdict,new_assembly[4])
                if la is None:
                    rejected_new_assemblies.append(new_assembly)
                    encountered_leaf_subsets.add(frozenset({frozenset({frozenset(n.leaf_names) for n in T(c).get_nonterminals()}) for c in new_assembly[0]}))
                elif len(accepted_assemblies) == num_requested_trees and la < accepted_assemblies[-1][4]:
                    rejected_new_assemblies.append(new_assembly)
                    encountered_leaf_subsets.add(frozenset({frozenset({frozenset(n.leaf_names) for n in T(c).get_nonterminals()}) for c in new_assembly[0]}))
                else:
                    leafsets = frozenset({frozenset({frozenset(n.leaf_names) for n in T(c).get_nonterminals()}) for c in new_assembly[0]})
                    if leafsets in encountered_leaf_subsets:
                        rejected_new_assemblies.append(new_assembly)
                    else:
                        encountered_leaf_subsets.add(leafsets)
    #                 for other_assembly in assemblies:
#                     if new_assembly[6] == other_assembly[6]:
#                         assert set(new_assembly[2]) == set(other_assembly[2])
#                         rejected_new_assemblies.append(new_assembly)
#                         break
#                 else:
#                     if new_assembly[6] in encountered_leaf_subsets:
#                         rejected_new_assemblies.append(new_assembly)
        new_assemblies = filter(lambda x: x not in rejected_new_assemblies,new_assemblies)
        new_this_iteration.extend(new_assemblies)
        assemblies.extend(new_assemblies)
        if still_in_play and len(assembly[0]) == 1 and not assembly[1]:
            accepted_this_iteration.append(assembly)
            accepted_assemblies.append(assembly)
            accepted_assemblies = sorted(accepted_assemblies,key=lambda x: x[4],reverse=True)
            if len(accepted_assemblies) > num_requested_trees:
                assert len(accepted_assemblies)-1 == num_requested_trees
                if accepted_assemblies[-1] in assemblies and accepted_assemblies[-1] not in discarded_this_iteration['cs']:
                    discarded_this_iteration['cs'].append(assemblies.pop(assemblies.index(accepted_assemblies.pop())))
                elif accepted_assemblies[-1] not in discarded_this_iteration['cs']:
                    discarded_this_iteration['cs'].append(accepted_assemblies.pop())
        for new_assembly in new_assemblies:
            if len(new_assembly[0]) == 1 and not new_assembly[1]:
                accepted_this_iteration.append(new_assembly)
                accepted_assemblies.append(new_assembly)
                accepted_assemblies = sorted(accepted_assemblies,key=lambda x: x[4],reverse=True)
                if len(accepted_assemblies) > num_requested_trees:
                    assert len(accepted_assemblies)-1 == num_requested_trees
                    if accepted_assemblies[-1] in assemblies and accepted_assemblies[-1] not in discarded_this_iteration['cs']:
                        discarded_this_iteration['cs'].append(assemblies.pop(assemblies.index(accepted_assemblies.pop())))
                    elif accepted_assemblies[-1] not in discarded_this_iteration['cs']:
                        discarded_this_iteration['cs'].append(accepted_assemblies.pop())
    print >>stderr,' '*110+'\r',
    if new_this_iteration:
        print >>stderr, '\t',len(new_this_iteration),"were added on this iteration"
    if discarded_this_iteration['u']:
        print >>stderr, '\t',len(discarded_this_iteration['u']),"were discarded because they became non-unique"
    if discarded_this_iteration['e']:
        print >>stderr, '\t',len(discarded_this_iteration['e']),"were discarded because there was no way to extend them"
    if discarded_this_iteration['cs']:
        print >>stderr, '\t',len(discarded_this_iteration['cs']),"were discarded because their current score is worse than the worst accepted score"
    if discarded_this_iteration['ps']:
        print >>stderr, '\t',len(discarded_this_iteration['ps']),"were discarded because their best possible final score is worse than the worst accepted score"
    if accepted_this_iteration:
        print >>stderr, '\t',len(accepted_this_iteration),"were completed and accepted on this iteration"
    if accepted_assemblies:
        print >>stderr, '\t'+"Best accepted score is",accepted_assemblies[0][4],"and worst accepted score is",accepted_assemblies[-1][4]
    print >>stderr,'\t',len(encountered_leaf_subsets),"unique leaf sets have been encountered"
    return accepted_assemblies


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
        self.unverified.pop(pair.leaves) # and is "verified", so pop it from unverified
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
    
    try:
      current_clades = self._nested_set_reprs
    except AttributeError:
      self._nested_set_reprs = [frozenset({c.nested_set_repr(),'r'})
                                for c in self.built_clades]
      current_clades = self._nested_set_reprs
    
    old_clades = frozenset(c for i,c in enumerate(current_clades)
                           if i not in indeces_to_skip)
    return frozenset({new_clade,old_clades})
  
  def best_case_with_extension(self,extension):
    if not hasattr(self,'_distances_to_root'):
      self._distances_to_root = {leaf:clade.trace_dist(leaf) for clade in self.built_clades
                                 for leaf in clade.leaf_names}
    if not hasattr(self,'_pairs_accounted_for'):
      self._pairs_accounted_for = {frozenset(pair) for clade in self.built_clades
                                   for pair in itertools.combinations(clade.leaf_names,2)}
    try:
      updated_distances_to_root = {leaf:(self._distances_to_root[leaf]+1 if leaf in
                                         self._distances_to_root else 1)
                               for pair in extension.consistent for leaf in pair}
      for leaf in self._distances_to_root:
        if leaf not in updated_distances_to_root:
          updated_distances_to_root[leaf] = self._distances_to_root[leaf]
      
      updated_pairs_accounted_for = {pair for pair in
                                     itertools.chain(self._pairs_accounted_for,
                                                     extension.consistent.iterkeys())}
    except AttributeError:
      updated_distances_to_root = dict(self._distances_to_root.iteritems())
      for leaf in extension.leaves:
        updated_distances_to_root[leaf] = 1
      updated_pairs_accounted_for = {pair for pair in self._pairs_accounted_for}
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
            if extension[1].score + self.score < min_score:
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
        build_in.score += math.log(pair.freq)
        updated_assemblies.append(build_in)
    return updated_assemblies
  
  def generate_extensions(self,encountered_assemblies,min_score=None):
    new_pairs,joins,attachments = self.find_extensions(encountered_assemblies,min_score=None)
    if any((new_pairs,joins,attachments)):
      # Return built extensions and any additional information
      return self.build_extensions(new_pairs, joins, attachments)
    else:
      # Somehow inform caller that this assembly is unextendable - it is the caller's responsibility to
      # remove this assembly from the container of active assemblies
      return None

def assemble_histtrees(pwhist,leaves_to_assemble,num_requested_trees=1000,freq_cutoff=0.9,max_iter=100000,processing_bite_size=10000):
    pwindiv = [(pair[0],score) for pair in [(f,[dd for i,dd in enumerate(ds)
                if sum(ddd[1] for ddd in ds[:i]) < freq_cutoff]) for f,ds in pwhist] for score in pair[1]]
    pairs_master = sorted(pwindiv,key=lambda x: (x[1][0],1-x[1][1]))
    assemblies = [[[],leaves_to_assemble,range(len(pairs_master)),[],0.0]]
    pwhist_as_dicts = {k:dict(v) for k,v in pwhist}
    iterations = 0
    accepted_assemblies = []
    encountered_leaf_subsets = set()
    enough_assemblies_accepted = False
    while any([len(asb[0]) > 1 or asb[1] for asb in assemblies]) and iterations < max_iter:
        new_accepted_assemblies = assembly_iteration(assemblies,pairs_master,pwhist_as_dicts,accepted_assemblies,
                                                 encountered_leaf_subsets,num_requested_trees,processing_bite_size,iterations)
        if new_accepted_assemblies is not accepted_assemblies:
            accepted_assemblies = new_accepted_assemblies
            with open("built_trees.txt",'w') as fh:
                for tree in accepted_assemblies:
                    treestr = StringIO()
                    T(tree[0][0]).write(treestr,'newick')
                    fh.write(str(tree[4])+'\t'+treestr.getvalue())
            print >>stderr, "Accepted assemblies have been updated"
        if not enough_assemblies_accepted and len(accepted_assemblies) >= num_requested_trees:
            assert len(accepted_assemblies) == num_requested_trees
            prev_assemblies_len = len(assemblies)
            assemblies = filter(lambda x: x[4] > accepted_assemblies[-1][4],assemblies)
            print >>stderr,"*** Reached",num_requested_trees,"accepted trees, dumping",prev_assemblies_len-len(assemblies),\
                            "assemblies with scores worse than the worst accepted tree ***"
            enough_assemblies_accepted = True
#         with open("trees_on_iter_"+str(iterations).zfill(2)+'.txt','w') as fh:
#             for asb in sorted(assemblies,key=lambda x: (len(x[1]),-x[4])):
#                 nwkclades = []
#                 for c in asb[0]:
#                     cladestr = StringIO()
#                     T(c).write(cladestr,'newick')
#                     nwkclades.append(cladestr.getvalue().strip().replace(':0.00000',''))
#                 fh.write('>>>'+str(asb[4])+'>>>'+'  '.join(nwkclades)+'\n')
#             print >>stderr,"Wrote current assemblies to file","trees_on_iter_"+str(iterations).zfill(2)+'.txt'
        iterations += 1
    return assemblies,accepted_assemblies
