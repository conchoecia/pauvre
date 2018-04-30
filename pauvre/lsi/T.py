# Binary search tree that holds status of sweep line. Only leaves hold values.
# Operations for finding left and right neighbors of a query point p and finding which segments contain p.
# Author: Sam Lichtenberg
# Email: splichte@princeton.edu
# Date: 09/02/2013

from pauvre.lsi.helper import *

ev = 0.00000001

class T:
	def __init__(self):
		self.root = Node(None, None, None, None)
	def contain_p(self, p):
		if self.root.value is None:
			return [[], []]
		lists = [[], []]
		self.root.contain_p(p, lists) 
		return (lists[0], lists[1])
	def get_left_neighbor(self, p):
		if self.root.value is None:
			return None
		return self.root.get_left_neighbor(p)
	def get_right_neighbor(self, p):
		if self.root.value is None:
			return None
		return self.root.get_right_neighbor(p)
	def insert(self, key, s):
		if self.root.value is None:
			self.root.left = Node(s, None, None, self.root)
			self.root.value = s 
			self.root.m = get_slope(s)
		else:
			(node, path) = self.root.find_insert_pt(key, s)
			if path == 'r':
				node.right = Node(s, None, None, node)
				node.right.adjust()
			elif path == 'l':
				node.left = Node(s, None, None, node)
			else:
				# this means matching Node was a leaf
				# need to make a new internal Node
				if node.compare_to_key(key) < 0 or (node.compare_to_key(key)==0 and node.compare_lower(key, s) < 1):
					new_internal = Node(s, None, node, node.parent)
					new_leaf = Node(s, None, None, new_internal)
					new_internal.left = new_leaf
					if node is node.parent.left:
						node.parent.left = new_internal
						node.adjust()
					else:
						node.parent.right = new_internal
				else:
					new_internal = Node(node.value, node, None, node.parent)
					new_leaf = Node(s, None, None, new_internal)
					new_internal.right = new_leaf
					if node is node.parent.left:
						node.parent.left = new_internal
						new_leaf.adjust()
					else:
						node.parent.right = new_internal
				node.parent = new_internal

	def delete(self, p, s):
		key = p
		node = self.root.find_delete_pt(key, s)
		val = node.value
		if node is node.parent.left:
			parent = node.parent.parent
			if parent is None:
				if self.root.right is not None:
					if self.root.right.left or self.root.right.right:
						self.root = self.root.right
						self.root.parent = None
					else:
						self.root.left = self.root.right
						self.root.value = self.root.right.value
						self.root.m = self.root.right.m
						self.root.right = None
				else:
					self.root.left = None
					self.root.value = None
			elif node.parent is parent.left:
				parent.left = node.parent.right
				node.parent.right.parent = parent
			else:
				parent.right = node.parent.right
				node.parent.right.parent = parent
		else:
			parent = node.parent.parent
			if parent is None:
				if self.root.left:
					# switch properties
					if self.root.left.right or self.root.left.left:
						self.root = self.root.left
						self.root.parent = None
					else:
						self.root.right = None
				else:
					self.root.right = None
					self.root.value = None
			elif node.parent is parent.left:
				parent.left = node.parent.left
				node.parent.left.parent = parent
				farright = node.parent.left
				while farright.right is not None:
					farright = farright.right
				farright.adjust()
			else:
				parent.right = node.parent.left
				node.parent.left.parent = parent
				farright = node.parent.left
				while farright.right is not None:
					farright = farright.right
				farright.adjust()
		return val

	def print_tree(self):
		self.root.print_tree()
class Node:
	def __init__(self, value, left, right, parent):
		self.value = value # associated line segment 
		self.left = left
		self.right = right
		self.parent = parent
		self.m = None
		if value is not None:
			self.m = get_slope(value)

	# compares line segment at y-val of p to p 
	# TODO: remove this and replace with get_x_at
	def compare_to_key(self, p):
		x0 = self.value[0][0]
		y0 = self.value[0][1]
		y1 = p[1]
		if self.m != 0 and self.m is not None:
			x1 = x0 - float(y0-y1)/self.m
			return compare_by_x(p, (x1, y1))
		else:
			x1 = p[0]
			return 0 
	
	def get_left_neighbor(self, p):
		neighbor = None
		n = self
		if n.left is None and n.right is None:
			return neighbor
		last_right = None
		found = False
		while not found:
			c = n.compare_to_key(p)
			if c < 1 and n.left:
				n = n.left
			elif c==1 and n.right:
				n = n.right
				last_right = n.parent
			else:
				found = True
		c = n.compare_to_key(p)
		if c==0:
			if n is n.parent.right:
				return n.parent
			else:
				goright = None
				if last_right:
					goright =last_right.left
				return self.get_lr(None, goright)[0]
		# n stores the highest-value in the left subtree
		if c==-1:
			goright = None
			if last_right:
				goright = last_right.left
			return self.get_lr(None, goright)[0]
		if c==1:
			neighbor = n
		return neighbor

	def get_right_neighbor(self, p):
		neighbor = None
		n = self
		if n.left is None and n.right is None:
			return neighbor
		last_left = None
		found = False
		while not found:
			c = n.compare_to_key(p)
			if c==0 and n.right:
				n = n.right
			elif c < 0 and n.left:
				n = n.left
				last_left = n.parent
			elif c==1 and n.right:
				n = n.right
			else:
				found = True
		c = n.compare_to_key(p)
		# can be c==0 and n.left if at root node
		if c==0:
			if n.parent is None:
				return None
			if n is n.parent.right:
				goleft = None
				if last_left:
					goleft = last_left.right
				return self.get_lr(goleft, None)[1]
			else:
				return self.get_lr(n.parent.right, None)[1]
		if c==1:
			goleft = None
			if last_left:
				goleft = last_left.right
			return self.get_lr(goleft, None)[1]
		if c==-1:
			return n
		return neighbor

	# travels down a single direction to get neighbors
	def get_lr(self, left, right):
		lr = [None, None]
		if left:
			while left.left:
				left = left.left
			lr[1] = left
		if right:
			while right.right:
				right = right.right
			lr[0] = right
		return lr
	
	def contain_p(self, p, lists):
		c = self.compare_to_key(p)
		if c==0:
			if self.left is None and self.right is None:
				if compare_by_x(p, self.value[1])==0:
					lists[1].append(self.value)
				else:
					lists[0].append(self.value)	
			if self.left:
				self.left.contain_p(p, lists)
			if self.right:
				self.right.contain_p(p, lists)
		elif c < 0:
			if self.left:
				self.left.contain_p(p, lists)
		else:
			if self.right:
				self.right.contain_p(p, lists)

	def find_insert_pt(self, key, seg):
		if self.left and self.right:
			if self.compare_to_key(key) == 0 and self.compare_lower(key, seg)==1:
				return self.right.find_insert_pt(key, seg)
			elif self.compare_to_key(key) < 1:
				return self.left.find_insert_pt(key, seg)
			else:
				return self.right.find_insert_pt(key, seg)	
		# this case only happens at root
		elif self.left:
			if self.compare_to_key(key) == 0 and self.compare_lower(key, seg)==1:
				return (self, 'r')
			elif self.compare_to_key(key) < 1:
				return self.left.find_insert_pt(key, seg)
			else:
				return (self, 'r')
		else:
			return (self, 'n')

	# adjusts stored segments in inner nodes
	def adjust(self):
		value = self.value
		m = self.m
		parent = self.parent
		node = self
		# go up left as much as possible
		while parent and node is parent.right:
			node = parent
			parent = node.parent
		# parent to adjust will be on the immediate right
		if parent and node is parent.left:
			parent.value = value
			parent.m = m
	
	def compare_lower(self, p, s2):
		y = p[1] - 10
		key = get_x_at(s2, (p[0], y))
		return self.compare_to_key(key)

	# returns matching leaf node, or None if no match
	# when deleting, you don't delete below--you delete above! so compare lower = -1.
	def find_delete_pt(self, key, value):
		if self.left and self.right:
			# if equal at this pt, and this node's value is less than the seg's slightly above this pt
			if self.compare_to_key(key) == 0 and self.compare_lower(key, value)==-1:
				return self.right.find_delete_pt(key, value)
			if self.compare_to_key(key) < 1:
				return self.left.find_delete_pt(key, value)
			else:
				return self.right.find_delete_pt(key, value)
		elif self.left:
			if self.compare_to_key(key) < 1:
				return self.left.find_delete_pt(key, value)
			else:
				return None
		# is leaf
		else:
			if self.compare_to_key(key)==0 and segs_equal(self.value, value):
				return self
			else:
				return None

	# also prints depth of each node
	def print_tree(self, l=0):
		l += 1
		if self.left:
			self.left.print_tree(l)
		if self.left or self.right:
			print('INTERNAL: {0}'.format(l))
		else:
			print('LEAF: {0}'.format(l))
		print(self)
		print(self.value)
		if self.right:
			self.right.print_tree(l)
