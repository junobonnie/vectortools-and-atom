# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 01:00:19 2021

@author: junob
"""
import math as m

class Tensor:
    def __init__(self, xx, xy, yx, yy):
        self.xx = xx
        self.xy = xy
        self.yx = yx
        self.yy = yy
        
    def __abs__(self):
        return self.xx*self.yy - self.xy*self.yx
    
    def __add__(self, other_Tensor):
        return Tensor(self.xx + other_Tensor.xx, self.xy + other_Tensor.xy,
                      self.yx + other_Tensor.yx, self.yy + other_Tensor.yy)
    
    def __sub__(self, other_Tensor):
        return Tensor(self.xx - other_Tensor.xx, self.xy - other_Tensor.xy,
                      self.yx - other_Tensor.yx, self.yy - other_Tensor.yy)
    
    def __mul__(self, num):
        return Tensor(self.xx*num, self.xy*num, self.yx*num, self.yy*num)
    
    def __rmul__(self, num):
        return Tensor(self.xx*num, self.xy*num, self.yx*num, self.yy*num)
    
    def __truediv__(self, num):
        return Tensor(self.xx/num, self.xy/num, self.yx/num, self.yy/num)
    
    def __str__(self):
        return 'Tensor(' + str(self.xx) + ', ' + str(self.xy) + ', ' + str(self.yx) + ', ' + str(self.yy) + ')'
    
    def dot(self, other_object):
        if isinstance(other_object, Tensor):
            return Tensor(self.xx*other_object.xx + self.xy*other_object.yx, self.xx*other_object.xy + self.xy*other_object.yy, 
                          self.yx*other_object.xx + self.yy*other_object.yx, self.yx*other_object.xy + self.yy*other_object.yy)
            
        elif isinstance(other_object, Vector):
            return Vector(self.xx*other_object.x + self.xy*other_object.y, self.yx*other_object.x + self.yy*other_object.y)
        
class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
    def __abs__(self):
        return m.sqrt(self.x**2 + self.y**2)
    
    def __add__(self, other_Vector):
        return Vector(self.x + other_Vector.x, self.y + other_Vector.y)
    
    def __sub__(self, other_Vector):
        return Vector(self.x - other_Vector.x, self.y - other_Vector.y)
    
    def __mul__(self, num):
        return Vector(self.x*num, self.y*num)
    
    def __rmul__(self, num):
        return Vector(self.x*num, self.y*num)
    
    def __truediv__(self, num):
        return Vector(self.x/num, self.y/num)
    
    def __str__(self):
        return 'Vector(' + str(self.x) + ', ' + str(self.y) + ')'
    
    def dot(self, other_Vector):
        return self.x*other_Vector.x + self.y*other_Vector.y
    
    def cross(self, other_Vector):
        return self.x*other_Vector.y - self.y*other_Vector.x
    
def SO2(theta):
    return Tensor(m.cos(theta), m.sin(theta), -m.sin(theta), m.cos(theta))