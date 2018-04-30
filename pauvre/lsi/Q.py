# Binary search tree that holds status of sweep line. Only leaves hold values.
# Operations for finding left and right neighbors of a query point p and finding which segments contain p.
# Author: Sam Lichtenberg
# Email: splichte@princeton.edu
# Date: 09/02/2013

from pauvre.lsi.helper import *

ev = 0.00000001

class Q:
	def __init__(self, key, value):
		self.key = key
		self.value = value
		self.left = None
		self.right = None
	
	def find(self, key):
		if self.key is None:
			return False
		c = compare_by_y(key, self.key)
		if c==0:
			return True
		elif c==-1:
			if self.left:
				self.left.find(key)
			else:
				return False
		else:
			if self.right:
				self.right.find(key)
			else:
				return False
	def insert(self, key, value):
		if self.key is None:
			self.key = key
			self.value = value
		c = compare_by_y(key, self.key)
		if c==0:
			self.value += value
		elif c==-1:
			if self.left is None:
				self.left = Q(key, value)
			else:
				self.left.insert(key, value)
		else:
			if self.right is None:
				self.right = Q(key, value)
			else:
				self.right.insert(key, value)
	# must return key AND value
	def get_and_del_min(self, parent=None):
		if self.left is not None:
			return self.left.get_and_del_min(self)
		else:
			k = self.key
			v = self.value
			if parent:
				parent.left = self.right
			# i.e. is root node
			else:
				if self.right:
					self.key = self.right.key
					self.value = self.right.value
					self.left = self.right.left
					self.right = self.right.right
				else:
					self.key = None
			return k,v
	
	def print_tree(self):
		if self.left:
			self.left.print_tree()
		print(self.key)
		print(self.value)
		if self.right:
			self.right.print_tree()
