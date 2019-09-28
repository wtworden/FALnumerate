
import snappy
import itertools
import operator
import csv
from os.path import expanduser


#### some of the below code relies on the plantri package, which is not installed in Sage by default. 


#### TOP-LEVEL CODE FOR GENERATING FAL CENSUSES TO FILE ----------------------------------------------------------------------------------------#####

## Note: HTP = half-twist partners
def FALnumerate(num_circles,return_DT_Codes=True,filename=None, identify=False,prime_only=False, group_into_HTP_classes=False, one_per_iso_class=False):
    DT_Codes = DTC_enumerate(num_circles,prime_only,group_into_HTP_classes)
    if not group_into_HTP_classes:
        DT_Codes = [DT_Codes]
    else:
        DT_Codes = list(DT_Codes)
    if filename:
        if group_into_HTP_classes:
            with open(filename,'w') as f:
                f.write('census of hyperbolic FALs with {} crossing circles, grouped by HTP class then homeomorphism type (non-rigorously).\n\n'.format(num_circles))
        else:
            with open(filename,'w') as f:
                f.write('census of hyperbolic FALs with {} crossing circles, grouped by homeomorphism type (non-rigorously).\n\n'.format(num_circles))
                f.write('total number of homeomorphism classes: {}\n'.format(len(isomorphism_classes)))
                f.write('total number of FAL DT codes (these are not necessarily distinct as links): {}\n\n'.format(len([0 for c in isomorphism_classes for M in c])))
    
    iso_class_groups = []
    for code_group in DT_Codes:
        classes = num_cusps_classes(code_group)
        exterior_classes = []
        for c in classes:
            ext_class = []
            for DT in c:
                M = snappy.DTcodec(list(DT)).link().exterior()
                ext_class.append(M)
            exterior_classes.append(ext_class)
        classes = []
        for e in exterior_classes:
            vol_classes = volume_classes(e)
            for vc in vol_classes:
                classes.append(vc)
        isomorphism_classes = []
        for c in classes:
            isomorphism_classes += iso_classes(c)
        iso_class_groups.append(isomorphism_classes)
        if identify and not group_into_HTP_classes:
            in_census = {}
            for c in isomorphism_classes:
                ident = c[0].identify()
                if ident:
                    names = tuple(M.name() for M in ident)
                    in_census[c[0]]=names
    
        if filename:
            with open(filename,'a') as f:
                f.write('new HTP class:\n')
                if identify and not group_into_HTP_classes:
                    f.write('FALs in SnapPy census:\n')
                    for n in in_census.values():
                        f.write('{}\n'.format(n))
                    f.write('\n')
    
                for c in isomorphism_classes:
                    f.write('Volume: {}\n'.format(c[0].volume()))
                    if one_per_iso_class:
                        f.write('{}\n'.format(c[0].DT_code()))
                    else:
                        for M in c:
                            f.write('{}\n'.format(M.DT_code()))
                f.write('\n\n')
    if group_into_HTP_classes:
        return None
    else:
        if return_DT_Codes:
            DT_classes = []
            for c in isomorphism_classes:
                DT_class = [M.DT_code() for M in c]
                DT_classes.append(DT_class)
            return DT_classes
        return isomorphism_classes, in_census

### generate a census of all FAL complements have n crossing circles, for all n in the list num_circles. Output is a spreadsheet which also has other 
### geometric and topological information.
def FALsheet(num_circles, path, filename):
    home=expanduser("~")
    with open('{}/Dropbox/research/active/Hidden_symmetries/computation/{}'.format(home,filename),'w') as csv_file:
        fieldnames = ['num crossing circles', 'num planar comps', 'volume', 'is_two_bridge', 'symmetry group of ext', 'homeomorphism class', 'census name','DT code']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
    for n in num_circles:
        n_iso_class = 0
        iso_classes, in_census = FALnumerate(n, return_DT_Codes=False, filename=None, identify=True, prime_only=False, group_into_HTP_classes=False, one_per_iso_class=False)
        for c in iso_classes:
            for ext in c:
                with open('{}/Dropbox/research/active/Hidden_symmetries/computation/{}'.format(home,filename),'a') as csv_file:
                    fieldnames = ['num crossing circles','num planar comps','volume','is_two_bridge','symmetry group of ext', 'homeomorphism class', 'census name','DT code']
                    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                    if c[0] in in_census:
                        names = in_census[c[0]]
                    else:
                        names = ''
                    iso_class = '{}-{}'.format(n,n_iso_class)
                    if ext.num_cusps()-n == 1:
                        N = ext.copy()
                        filling = [(0,0)]+[(1,1) for i in range(n)]
                        N.dehn_fill(filling)
                        M = N.filled_triangulation()
                        is_two_bridge = M.is_two_bridge()
                    else:
                        is_two_bridge = ''
                    result = {'num crossing circles':n,'num planar comps':ext.num_cusps()-n,'volume':ext.volume(), 'is_two_bridge':is_two_bridge,'symmetry group of ext':ext.symmetry_group(), 'homeomorphism class':iso_class, 'census name':names,'DT code':ext.DT_code()}
                    writer.writerow(result)
            n_iso_class += 1

#### ----------------------------------------------------------------------------------------------------------------------------------------------#####




#### HELPERS FOR ENUMERATION FUNCTIONS ------------------------------------------------------------------------------------------------------------#####


### partition a set of link exteriors into subsets, such that all members of a subset have the same (approximate) volume.
def volume_classes(link_exteriors):
    classes = []
    exts_w_vol = [(l,int(l.volume()*1000000)) for l in link_exteriors]
    exts_w_vol.sort(key=lambda x: x[1])
    for key,group in itertools.groupby(exts_w_vol,operator.itemgetter(1)):
        L_class = list(group)
        c = [L[0] for L in L_class]
        classes.append(c)
    return classes

### partition a set of link given by DT codes into subsets, such that all members of a subset have the same number of cusps.
def num_cusps_classes(DT_Codes):
    classes = []
    DTs_w_lengths = [(C,len(C)) for C in DT_Codes]
    DTs_w_lengths.sort(key=lambda x: x[1])
    for key,group in itertools.groupby(DTs_w_lengths,operator.itemgetter(1)):
        DT_class = list(group)
        c = [DT[0] for DT in DT_class]
        classes.append(c)
    return classes

### partition a set of link exteriors into isomorphism classes.
def iso_classes(link_exteriors):
    classes = []
    c = []
    c.append(link_exteriors.pop())
    classes.append(c)
    while link_exteriors:
        L = link_exteriors[-1]
        for c in classes:
            if L.is_isometric_to(c[0]):
                c.append(link_exteriors.pop())
                L = None
                break
        if L:
            classes.append([link_exteriors.pop()])
    return classes

#### ----------------------------------------------------------------------------------------------------------------------------------####





#### CODE FOR COMPUTING DT CODES ------------------------------------------------------------------------------------------------------####

## Square is a class for storing the crossing labels for DT code computation. For each marked edge  of the dual nerve N we get a square. 
## Marked edges of the dual nerve are considered to be oriented to point toward the vertex of higher index, and the square associated
## to such an edge is oriented so that the edge point from bottom to top of the square. For each square we have
## four crossings coming from the crossing circle, called upper_right, upper_left, lower_right, and lower_left, and one additional crossing
## if there is a half-twist, called middle. Each of these has associated to it a tuple, which at initiation is [0,0]. To compute the DT code 
## we walk around the link and label each crossing as we enter it, once for each time. These labels are stored in the tuple. Each time a 0
## in one of these tuples changes to the crossing label (which is non-zero), we reduce self.num_zeros by one. When self.num_zeros==0, we 
## know that all crossing associated to this square have been visited twice.

#      ____________________________________________________o                                ____________________________________________________o 
#      |    \                                          /   |                                |    \                                          /   |  
#      |      \                                      /     |                                |      \                                      /     |  
#      |        \                                  /       |                                |        \                                  /       |  
#      |          \                              /         |      N:                        |          \                              /         |  
#      |            \                          /           |                                |            \                          /           |  
#      |              \                      /             |         \         /            |              \                      /             |  
#      |            UL  \                  / UR            |           \     /              |            UL  \                  / UR            |  
#      |              _______________________              |             \ /                |              _______________________              |  
#      |            /      \           /      \            |              | j               |            /      \           /      \            |  
#      |          /          \       /          \          |              |                 |          /          \       /          \          |  
#      |        /              \   /              \        |              ^                 |        /             |     |             \        |  
#      |       |                 /   M              |      |  <----'X'    |      'I'---->   |       |              |     |              |       | 
#      |        \              /   \              /        |              | e               |        \             |     |             /        |  
#      |          \      LL  /       \  LR      /          |              |                 |          \      LL  /       \  LR      /          |  
#      |            \______ /_________\ ______/            |              |                 |            \______ /_________\ ______/            |   
#      |                  /             \                  |              |                 |                  /             \                  |  
#      |                /                 \                |              | i               |                /                 \                |  
#      |              /                     \              |             / \                |              /                     \              |  
#      |            /                         \            |           /     \              |            /                         \            |  
#      |          /                             \          |         /         \            |          /                             \          |  
#      |        /                                 \        |                                |        /                                 \        |  
#      |      /                                     \      |                                |      /                                     \      |  
#      |    /                                         \    |                                |    /                                         \    |  
#      |  /                                             \  |                                |  /                                             \  |  
#      |___________________________________________________|                                |___________________________________________________|





class Square:
    def __init__(self,i,j,label):
        self.index = None
        self.edge = (i,j)
        self.upper_right = [0,0]
        self.upper_left = [0,0]
        self.lower_right = [0,0]
        self.lower_left = [0,0]
        self.middle = None if label == 'I' else [0,0]
        self.non_empty = 0 #empty
        self.label = label
        self.full = 0
        self.num_zeros = 8 if label == 'I' else 10

    def corners(self, corner):
        if corner == 'upper_left':
            return self.upper_left
        elif corner == 'upper_right':
            return self.upper_right
        elif corner == 'lower_left':
            return self.lower_left
        elif corner == 'lower_right':
            return self.lower_right
        elif corner == 'middle':
            return self.middle



## A (dual) nerve corresponds to a link, so there is a corresponding DT code. This is useful because it can be plugged in to SnapPy. 
## For an explanation of how (in general) to compute DT codes, see http://katlas.org/wiki/DT_(Dowker-Thistlethwaite)_Codes.
## To compute the DT code we walk along the link, labelling crossings as we go as describe at the above link. labellings are stores in 
## Squares. If we know via which corner we enter a square, then the corner we leave be is determined by the labelling of the square by
## 'X' or 'I' (e.g., if we enter the upper-left corner of an 'X' labelled square, we'll leave via the lower-right corner). When leaving a 
## square, we must determine which square we will enter next, and via which corner. Since we know how the square is align with the edge, this
## is determined by the dual nerve N. For example, if we leave the upper right corner, then we look at the vertex v of e of larger index, and 
## follow the edge e0 meeting v that is adjacent to the face to the right of e. At the other end of this edge is a vertex v', and some marked edge 
## e' meets v'. The square associated to e' is the next square, and the corner we enter this square depends on whether v' is the vertex of e'
## of largest index, and whether it also bounds the face bounded by e and e0. This is all pretty understandable if you draw a picture of an FAL 
## overlaid with a picture of the dual nerve, and put a square on marked edges, enclosing the associated crossing circle.
def compute_DT_code(dual_nerve):
    N = dual_nerve
    E = N.get_embedding()
    squares = [Square(e[0],e[1],e[2]) for e in N.edges() if e[2] != None]
    for i in range(len(squares)):
        squares[i].index = i
    unfilled_squares = [i for i in range(len(squares))]
    components = []
    count = 1
    while unfilled_squares:
        S = sorted([s for s in squares if s.full == 0],key=lambda x: x.non_empty)[-1]
        if not S.non_empty:
            crossing_labels = []
            count, unfilled_squares, crossing_labels = chase_planar(squares, unfilled_squares, crossing_labels, E, S, 'upper_right', count)
            components.append(crossing_labels)            
            assert S.non_empty == 1
            assert S.full == 0
        else:
            if ( S.upper_right[0] and S.upper_right[0] % 2 == count % 2 ) or ( S.upper_left[0] and S.upper_left[0] % 2 != count % 2 ):
                crossing_labels = []

                S.lower_right[1] = count
                crossing_labels.append(count); count += 1; S.num_zeros -= 1
                
                S.lower_left[1] = count
                crossing_labels.append(count); count += 1; S.num_zeros -= 1
                
                S.upper_left[1] = count if ( count % 2 == 1 ) else -count
                crossing_labels.append(count); count += 1; S.num_zeros -= 1
                
                S.upper_right[1] = count if ( count % 2 == 1 ) else -count
                crossing_labels.append(count); count += 1; S.num_zeros -= 1
                
                components.append(crossing_labels)
            else:
                crossing_labels = []
                S.upper_right[1] = count if ( count % 2 == 1 ) else -count
                crossing_labels.append(count); count += 1; S.num_zeros -= 1

                S.lower_right[1] = count
                crossing_labels.append(count); count += 1; S.num_zeros -= 1
                
                S.lower_left[1] = count
                crossing_labels.append(count); count += 1; S.num_zeros -= 1
                
                S.upper_left[1] = count if ( count % 2 == 1 ) else -count
                crossing_labels.append(count); count += 1; S.num_zeros -= 1
                
                components.append(crossing_labels)
            if ( ( not S.upper_right[0] ) and  ( S.upper_right[1] % 2 != count % 2 ) ):
                crossing_labels = []
                count, unfilled_squares, crossing_labels = chase_planar(squares, unfilled_squares, crossing_labels, E, S,'upper_right', count)
                components.append(crossing_labels)
            elif ( ( not S.upper_right[0] ) and ( S.upper_right[1] % 2 == count % 2 ) ):
                crossing_labels = []
                S, corner = next_square(squares, unfilled_squares, E, S, 'upper_right')
                assert S.corners(corner)[0] == 0
                count, unfilled_squares, crossing_labels = chase_planar(squares, unfilled_squares, crossing_labels, E, S, corner, count)
                components.append(crossing_labels)
            elif ( ( not S.upper_left[0] ) and ( S.upper_left[1] % 2 != count % 2 ) ):
                crossing_labels = []
                count, unfilled_squares, crossing_labels = chase_planar(squares, unfilled_squares, crossing_labels, E, S,'upper_left', count)
                components.append(crossing_labels)
            elif ( ( not S.upper_left[0] ) and ( S.upper_left[1] % 2 == count % 2 ) ):
                S, corner = next_square(squares, unfilled_squares, E, S, 'upper_left')
                crossing_labels = []
                assert S.corners(corner)[0] == 0
                count, unfilled_squares, crossing_labels = chase_planar(squares, unfilled_squares, crossing_labels, E, S, corner, count)
                components.append(crossing_labels)
            elif S.num_zeros == 0:
                S.full = 1
                unfilled_squares.remove(S.index)
    tuples = []
    for s in squares:
        for corner in ['upper_left','upper_right','lower_left','lower_right']:
            l = s.corners(corner)
            if l[0] % 2 == 0:
                l.reverse()
            for i in range(len(components)):
                if l[0] in components[i]:
                    l.append(i)
                    break
            tuples.append(tuple(l))
        if s.label =='X':
            l = s.middle
            if l[0] % 2 == 0:
                l.reverse()
            for i in range(len(components)):
                if l[0] in components[i]:
                    l.append(i)
                    break
            tuples.append(tuple(l))
    tuples.sort(key=lambda x: x[0])
    DT_tuples = []
    while tuples:
        comp = tuples[0][2]
        DT_tup = [tuples[0][1]]
        tuples.remove(tuples[0])
        while tuples and tuples[0][2] == comp:
            DT_tup.append(tuples[0][1])
            tuples.remove(tuples[0])
        DT_tuples.append(tuple(DT_tup))

    return DT_tuples

def get_square(i, j, squares):
    i,j = sorted([i,j])
    for s in squares:
        if s.edge == (i,j):
            return s
    return 0

def is_square(i, j, squares):
    i,j = sorted([i,j])
    for s in squares:
        if s.edge == (i,j):
            return s
    return 0

## if we enter a square with the given label via the given corner, returns the corner through which we leave.
def next_corner(corner, label):
    if (corner, label) == ('upper_right','I'):
        return 'lower_right'
    elif (corner, label) == ('upper_right', 'X'):
        return 'lower_left'
    if (corner, label) == ('upper_left','I'):
        return 'lower_left'
    elif (corner, label) == ('upper_left', 'X'):
        return 'lower_right'
    if (corner, label) == ('lower_right','I'):
        return 'upper_right'
    elif (corner, label) == ('lower_right', 'X'):
        return 'upper_left'
    if (corner, label) == ('lower_left','I'):
        return 'upper_left'
    elif (corner, label) == ('lower_left', 'X'):
        return 'upper_right'

## find the next square and the corner through which we enter it.
def next_square(squares, unfilled_squares, embedding, square, corner):
    S = square
    E = embedding
    if corner in ['upper_right', 'upper_left']:
        pivot1 = S.edge[1]
        E_loc = E[pivot1]
        new_sq_vert1 = E_loc[E_loc.index(S.edge[0]) - 1] if corner == 'upper_right' else E_loc[(E_loc.index(S.edge[0]) + 1) % len(E_loc)]
        E_loc = E[new_sq_vert1]
        new_sq_vert2 = E_loc[E_loc.index(S.edge[1]) - 1]

    elif corner in ['lower_right', 'lower_left']:
        pivot1 = S.edge[0]
        E_loc = E[pivot1]
        new_sq_vert1 = E_loc[E_loc.index(S.edge[1]) - 1] if corner == 'lower_left' else E_loc[(E_loc.index(S.edge[1]) + 1) % len(E_loc)]
        E_loc = E[new_sq_vert1]
        new_sq_vert2 = E_loc[E_loc.index(S.edge[0]) - 1]
    S0 = get_square(new_sq_vert1, new_sq_vert2, squares)
    if S0:
        next_sq = get_square(new_sq_vert1, new_sq_vert2, squares)
        next_corner = 'upper_left' if new_sq_vert1 > new_sq_vert2 else 'lower_right'
    else:
        new_sq_vert2 = E_loc[(E_loc.index(pivot1) + 1) % len(E_loc)]
        next_sq = get_square(new_sq_vert1, new_sq_vert2, squares)                
        next_corner = 'upper_right' if new_sq_vert1 > new_sq_vert2 else 'lower_left'
    return next_sq, next_corner
       
## if we label an overcrossing by an ever number, it gets a negative sign (see http://katlas.org/wiki/DT_(Dowker-Thistlethwaite)_Codes).
def sign(corner,parity):
    if parity == 1:
        return 1
    elif corner in ['lower_left','lower_right']:
        return -1
    else:
        return 1


def chase_planar(squares, unfilled_squares, crossing_labels, embedding, square, corner, count):
    S = square
    E = embedding
    S.corners(corner)[0] = count*sign(corner,count % 2); crossing_labels.append(count); count += 1; S.num_zeros -= 1
    if S.label == 'X':
        if corner in ['upper_left','lower_right']:
            S.middle[0] = count; crossing_labels.append(count); count += 1; S.num_zeros -= 1
        elif corner in ['upper_right','lower_left']:
            mult = -1 if ( count % 2 == 0 ) else 1
            S.middle[1] = count*mult; crossing_labels.append(count); count += 1; S.num_zeros -= 1

    corner = next_corner(corner, S.label)
    S.corners(corner)[0] = count*sign(corner,count % 2); crossing_labels.append(count); count += 1; S.num_zeros -= 1
    S.non_empty = 1
    if S.num_zeros == 0:
        S.full = 1
        unfilled_squares.remove(S.index)
    S, corner = next_square(squares, unfilled_squares, embedding, S, corner)
    if S.corners(corner)[0] == 0:
        count, unfilled_squares, crossing_labels = chase_planar(squares, unfilled_squares, crossing_labels, embedding, S, corner, count)
    
    return count, unfilled_squares, crossing_labels


## Note: HTP = half-twist partners
def DTC_enumerate(num_circles,prime_only=False, group_into_HTP_classes=False):
    dual_nerves = marked_dual_nerves(num_circles,prime_only,group_into_HTP_classes)
    DT_Codes = set([])
    for i in range(len(dual_nerves)):
        N = dual_nerves[i]
        if group_into_HTP_classes:
            DT_Code = tuple(tuple(compute_DT_code(dual_nerve)) for dual_nerve in N)
        else:
            DT_Code = tuple(compute_DT_code(N))
        DT_Codes.add(DT_Code)
    return DT_Codes




### CODE FOR COMPUTING MARKED DUAL NERVES -----------------------------------------------------------------------------------------###






## determine if a marked dual nerve corresponds with a prime FAL. For a nerve, non-primeness is equivalent to having a 3-cycle that does not 
## bound a triangle. Thus for a dual nerve, non-primeness is equivalent to having a triple of faces that are pairwise adjacent via three edges
## that do not all meet a common vertex. Note that one of these three edges must be a marked edge, and none of these faces can be triangles.
def is_prime(marked_dual_nerve):
    N = marked_dual_nerve
    F = N.faces()
    F_sorted = [[tuple(sorted(e)) for e in face] for face in F]
    E_marked = [e for e in N.edges() if e[2] != None] ## get all edges that correspond to crossing circles
    for e in E_marked:  ## for each marked edge, see if the adjacent faces have another adjacency
        i,j = e[0],e[1]  ##get the vertices of e
        e_face0, e_face1 = [face for face in F_sorted if (i,j) in face] ## get the faces adjacent to e
        if len(e_face0) > 3 and len(e_face1) > 3:  ## we only need to go further if neither face is a triangle..
            disjoint_edges0 = [e for e in e_face0 if ( ( i not in e ) and ( j not in e ) and ( e not in E_marked) )] ## get edges of e_face0 that are disjoint from e, and unmarked
            disjoint_edges1 = [e for e in e_face1 if ( ( i not in e ) and ( j not in e ) and ( e not in E_marked) )]
            for face in F_sorted:  ## check if there is a face that is adjacent to both e_face0 and e_face1
                if ( len(set(face).intersection(set(disjoint_edges0)))>0 and len(set(face).intersection(set(disjoint_edges1)))>0 ):
                    return False
    return True


## Note: HTP = half-twist partners
## Generate all possible marked dual nerves for FALs having num_circles crossing circles. A marked dual nerve is a the dual graph to the nerve, with 
## edges corresponding to crossing circles marked according to whether there is a half-twist at the crossing circle. If group_into_HTP_classes is True, 
## then all of the resulting dual nerves which differ only by their markings will be grouped together. If prime_only is true, then only 
## dual nerves of prime FALs will be returned. A prime FAL is one that is not a belted sum of other FALs, see Morgan--Spyropoulos--Trapp--Ziegler, 
## "Belted Sum Decomposition of Fully Augmented Links".
## Returns a list (with sublists if group_into_HTP_classes=True), each element of which is a graph G, with edges corresponding to crossing circles labelled with an 'X' or an 'I'.
def marked_dual_nerves(num_circles, prime_only=False, group_into_HTP_classes=False):
    dual_nerves_tups = enumerate_dualnerves(num_circles)
    crossing_words = get_words(num_circles)
    dual_nerves_w_crossings = []
    for tup in dual_nerves_tups:
        G = tup[0]
        for n in tup[1]:
            HTP_class = []
            for w in crossing_words:
                G_copy = G.copy()
                letter = 0
                for i in range(2*num_circles):
                    if letter < num_circles:
                        for j in range(i,2*num_circles):
                            if n[i][j]==-1:
                                G_copy.set_edge_label(i+1,j+1,w[letter])
                                letter += 1
                    else:
                        if group_into_HTP_classes:
                            HTP_class.append(G_copy)
                        else:
                            dual_nerves_w_crossings.append(G_copy)
                        break
            if group_into_HTP_classes:
                dual_nerves_w_crossings.append(HTP_class)
    if prime_only:
        if group_into_HTP_classes:
            prime_dual_nerves = [[N for N in HTP_class if is_prime(N)] for HTP_class in dual_nerves_w_crossings]
        else:
            prime_dual_nerves = [N for N in dual_nerves_w_crossings if is_prime(N)]
        return prime_dual_nerves
    return dual_nerves_w_crossings


## Generate all words in the letters 'X' and 'I' of length length. These words are used to mark dual nerves with half-twist 
## information. An 'X' corresponds to a half-twist, while an 'I' corresponds to no twist.
def get_words(length):
    words = []
    for l in ['X','I']:
        word = l
        next_letter(word,length,words)
    return words

def next_letter(word,length,words):
    if len(word)<length:
        for l in ['X','I']:
            new_word = ''.join([let for let in word])
            new_word += l
            next_letter(new_word,length,words)
    else:
        words.append(word)

## Enumerate all dual nerves of FALs with num_circles crossing circles. A nerve is a triangulation of S^2, and a dual nerve is
## the dual of such a triangulation, which is a trivalent graph. We use Sage function graphs.triangulations(), which uses the 
## plantri package. This is an optional package, so it has to be installed in Sage. 
## The edges of dual nerves are labelled by label_edges(), which for each nerves, finds every way of labelling the edges
## so that for each vertex, exactly two of the edges incident on it are labelled by 1, and one of the incident edges is labelled
## -1. A label of -1 indicated a crossing circle edge. Returns a list of tuples. Each tuple contains the underlying dual graph G,
## along with a matrix which is the same as the adjacency matrix of G, except that entries indicating an edge corresponding to
## a crossing circle are -1 instead of 1.
def enumerate_dualnerves(num_circles):
    c = num_circles
    nerve_num_verts = c+2
    dual_graphs = graphs.triangulations(nerve_num_verts,dual=True)
    dual_nerves_tups = []
    for G in dual_graphs:
        dual_nerves = []
        A = G.adjacency_matrix()

        B = [list(A[i][j] for j in range(A.dimensions()[1])) for i in range(A.dimensions()[0])]
        label_edges(B, dual_nerves)
        dual_nerves_tups.append((G,dual_nerves))
    return dual_nerves_tups


def label_edges(B, dual_nerves):
    unlabelled_rows = [i for i in range(len(B)) if sum(B[i])!=1]
    if unlabelled_rows:
        i = unlabelled_rows[0]
        good_ones = [k for k in range(len(B)) if B[i][k]==1 and k in unlabelled_rows]
        for j in good_ones:
            B_fork = [[B[k][l] for l in range(len(B))] for k in range(len(B))]
            B_fork[i][j] = -1
            B_fork[j][i] = -1
            label_edges(B_fork, dual_nerves)
    else:
        dual_nerves.append(B)      



