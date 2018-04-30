# Helper functions for use in the lsi implementation.

ev = 0.0000001
# floating-point comparison
def approx_equal(a, b, tol):
	return abs(a - b) < tol

# compares x-values of two pts
# used for ordering in T
def compare_by_x(k1, k2):
	if approx_equal(k1[0], k2[0], ev):
		return 0
	elif k1[0] < k2[0]:
		return -1
	else:
		return 1

# higher y value is "less"; if y value equal, lower x value is "less"
# used for ordering in Q
def compare_by_y(k1, k2):
	if approx_equal(k1[1], k2[1], ev):
		if approx_equal(k1[0], k2[0], ev):
			return 0
		elif k1[0] < k2[0]:
			return -1
		else:
			return 1
	elif k1[1] > k2[1]:
		return -1
	else:
		return 1

# tests if s0 and s1 represent the same segment (i.e. pts can be in 2 different orders)
def segs_equal(s0, s1):
	x00 = s0[0][0]
	y00 = s0[0][1]
	x01 = s0[1][0]
	y01 = s0[1][1]
	x10 = s1[0][0]
	y10 = s1[0][1]
	x11 = s1[1][0]
	y11 = s1[1][1]
	if (approx_equal(x00, x10, ev) and approx_equal(y00, y10, ev)):
		if (approx_equal(x01, x11, ev) and approx_equal(y01, y11, ev)):
			return True
	if (approx_equal(x00, x11, ev) and approx_equal(y00, y11, ev)):
		if (approx_equal(x01, x10, ev) and approx_equal(y01, y10, ev)):
			return True
	return False

# get m for a given seg in (p1, p2) form
def get_slope(s):
	x0 = s[0][0]
	y0 = s[0][1]
	x1 = s[1][0]
	y1 = s[1][1]
	if (x1-x0)==0:
		return None
	else:
		return float(y1-y0)/(x1-x0)

# given a point p, return the point on s that shares p's y-val
def get_x_at(s, p):
	m = get_slope(s)
	# TODO: this should check if p's x-val is octually on seg; we're assuming
	# for now that it would have been deleted already if not 
	if m == 0: # horizontal segment
		return p
	# ditto; should check if y-val on seg
	if m is None: # vertical segment
		return (s[0][0], p[1])
	x1 = s[0][0]-(s[0][1]-p[1])/m
	return (x1, p[1])

# returns the point at which two line segments intersect, or None if no intersection.
def intersect(seg1, seg2):
	p = seg1[0]
	r = (seg1[1][0]-seg1[0][0], seg1[1][1]-seg1[0][1])
	q = seg2[0]
	s = (seg2[1][0]-seg2[0][0], seg2[1][1]-seg2[0][1])
	denom = r[0]*s[1]-r[1]*s[0]
	if denom == 0:
		return None
	numer = float(q[0]-p[0])*s[1]-(q[1]-p[1])*s[0]
	t = numer/denom
	numer = float(q[0]-p[0])*r[1]-(q[1]-p[1])*r[0]
	u = numer/denom
	if (t < 0 or t > 1) or (u < 0 or u > 1):
		return None
	x = p[0]+t*r[0]
	y = p[1]+t*r[1]
	return (x, y)


