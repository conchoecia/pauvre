# Implementation of the Bentley-Ottmann algorithm, described in deBerg et al, ch. 2.
# See README for more information.
# Author: Sam Lichtenberg
# Email: splichte@princeton.edu
# Date: 09/02/2013

from pauvre.lsi.Q import Q
from pauvre.lsi.T import T
from pauvre.lsi.helper import *

# "close enough" for floating point
ev = 0.00000001

# how much lower to get the x of a segment, to determine which of a set of segments is the farthest right/left
lower_check = 100

# gets the point on a segment at a lower y value.
def getNextPoint(p, seg, y_lower):
	p1 = seg[0]
	p2 = seg[1]
	if (p1[0]-p2[0])==0:
		return (p[0]+10, p[1])
	slope = float(p1[1]-p2[1])/(p1[0]-p2[0])
	if slope==0:
		return (p1[0], p[1]-y_lower)
	y = p[1]-y_lower
	x = p1[0]-(p1[1]-y)/slope
	return (x, y)

"""
for each event point:
	U_p = segments that have p as an upper endpoint
	C_p = segments that contain p
	L_p = segments that have p as a lower endpoint
"""
def handle_event_point(p, segs, q, t, intersections):
	rightmost = (float("-inf"), 0)
	rightmost_seg = None
	leftmost = (float("inf"), 0) 
	leftmost_seg = None

	U_p = segs
	(C_p, L_p) = t.contain_p(p)
	merge_all = U_p+C_p+L_p
	if len(merge_all) > 1:
		intersections[p] = []
		for s in merge_all:
			intersections[p].append(s) 
	merge_CL = C_p+L_p
	merge_UC = U_p+C_p
	for s in merge_CL:
		# deletes at a point slightly above (to break ties) - where seg is located in tree
		# above intersection point
		t.delete(p, s)
	# put segments into T based on where they are at y-val just below p[1]
	for s in merge_UC:
		n = getNextPoint(p, s, lower_check) 
		if n[0] > rightmost[0]:
			rightmost = n 
			rightmost_seg = s
		if n[0] < leftmost[0]:
			leftmost = n
			leftmost_seg = s
		t.insert(p, s)

	# means only L_p -> check newly-neighbored segments
	if len(merge_UC) == 0:
		neighbors = (t.get_left_neighbor(p), t.get_right_neighbor(p))
		if neighbors[0] and neighbors[1]:
			find_new_event(neighbors[0].value, neighbors[1].value, p, q)
			
	# of newly inserted pts, find possible intersections to left and right
	else:
		left_neighbor = t.get_left_neighbor(p)
		if left_neighbor:
			find_new_event(left_neighbor.value, leftmost_seg, p, q)
		right_neighbor = t.get_right_neighbor(p)
		if right_neighbor:
			find_new_event(right_neighbor.value, rightmost_seg, p, q)

def find_new_event(s1, s2, p, q):
	i = intersect(s1, s2)
	if i:
		if compare_by_y(i, p) == 1:
			if not q.find(i):
				q.insert(i, [])
	
# segment is in ((x, y), (x, y)) form
# first pt in a segment should have higher y-val - this is handled in function
def intersection(S):
	s0 = S[0]
	if s0[1][1] > s0[0][1]:
		s0 = (s0[1], s0[0])
	q = Q(s0[0], [s0])
	q.insert(s0[1], [])
	intersections = {}
	for s in S[1:]:
		if s[1][1] > s[0][1]:
			s = (s[1], s[0])
		q.insert(s[0], [s])
		q.insert(s[1], [])
	t = T()
	while q.key:
		p, segs = q.get_and_del_min()
		handle_event_point(p, segs, q, t, intersections)
	return intersections

