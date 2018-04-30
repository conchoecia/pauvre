# Test file for lsi.
# Author: Sam Lichtenberg
# Email: splichte@princeton.edu
# Date: 09/02/2013

from lsi import intersection
import random
import time, sys
from helper import *

ev = 0.00000001

def scale(i):
	return float(i)

use_file = None
try:
	use_file = sys.argv[2]
except:
	pass

if not use_file:
	S = [] 
	for i in range(int(sys.argv[1])):
		p1 = (scale(random.randint(0, 1000)), scale(random.randint(0, 1000)))
		p2 = (scale(random.randint(0, 1000)), scale(random.randint(0, 1000)))
		s = (p1, p2)
		S.append(s)
	f = open('input', 'w')
	f.write(str(S))
	f.close()

else:
	f = open(sys.argv[2], 'r')
	S = eval(f.read())

intersections = []
seen = []
vs = False
hs = False
es = False
now = time.time()
for seg1 in S:
	if approx_equal(seg1[0][0], seg1[1][0], ev):
		print 'VERTICAL SEG'
		print ''
		print ''
		vs = True
	if approx_equal(seg1[0][1], seg1[1][1], ev):
		print 'HORIZONTAL SEG'
		print ''
		print ''
		hs = True
	for seg2 in S:
		if seg1 is not seg2 and segs_equal(seg1, seg2):
			print 'EQUAL SEGS'
			print ''
			print ''
			es = True
		if seg1 is not seg2 and (seg2, seg1) not in seen:
			i = intersect(seg1, seg2)
			if i:
				intersections.append((i, [seg1, seg2]))
		#		xpts = [seg1[0][0], seg1[1][0], seg2[0][0], seg2[1][0]]
		#		xpts = sorted(xpts)
		#		if (i[0] <= xpts[2] and i[0] >= xpts[1]:
		#			intersections.append((i, [seg1, seg2]))
				seen.append((seg1, seg2))
later = time.time()
n2time = later-now
print "Line sweep results:"
now = time.time()
lsinters = intersection(S)
inters = []
for k, v in lsinters.iteritems():
	#print '{0}: {1}'.format(k, v)
	inters.append(k)
#	inters.append(v)
later = time.time()
print 'TIME ELAPSED: {0}'.format(later-now)
print "N^2 comparison results:"
pts_seen = []
highestseen = 0
for i in intersections:
	seen_already = False
	seen = 0
	for p in pts_seen:
		if approx_equal(i[0][0], p[0], ev) and approx_equal(i[0][1], p[1], ev):
			seen += 1
			seen_already = True
	if seen > highestseen:
		highestseen = seen
	if not seen_already:
		pts_seen.append(i[0])
	in_k = False
	for k in inters:
		if approx_equal(k[0], i[0][0], ev) and approx_equal(k[1], i[0][1], ev):
			in_k = True
	if in_k == False:
		print 'Not in K: {0}: {1}'.format(i[0], i[1])
#	print i
print highestseen
print 'TIME ELAPSED: {0}'.format(n2time)
#print 'Missing from line sweep but in N^2:'
#for i in seen:
#	matched = False
print len(lsinters)
print len(pts_seen)
if len(lsinters) != len(pts_seen):
	print 'uh oh!'
