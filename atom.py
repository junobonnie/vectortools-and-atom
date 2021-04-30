# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 14:16:50 2021

@author: junob
"""
import math as m
from vectortools import *
import pygame as pg

class Wall:
    def __init__(self, width, height, theta, pos, color):
        self.width = width
        self.height = height
        self.theta = theta
        self.pos = pos
        self.color = color
        self.O = pos + SO2(theta).dot(Vector(width/2, height/2))
        
        self.A = pos
        self.B = pos + SO2(theta).dot(Vector(width, 0))
        self.C = pos + SO2(theta).dot(Vector(width, height))
        self.D = pos + SO2(theta).dot(Vector(0, height))
        
    def __str__(self):
        return ('Wall(width = ' + str(self.width) + ', height = ' + str(self.height) + 
                ', theta = ' + str(self.theta) + 
                ', pos(' + str(self.pos.x) + ', ' + str(self.pos.y) + ')' + 
                ', color = ' + str(self.color) + ')')
    
    def cross_3(self, A_Vector, B_Vector, C_Vector):
        AB = A_Vector - B_Vector
        AC = A_Vector - C_Vector
        if AB.cross(AC) > 0:
            return 1
        elif AB.cross(AC) < 0:
            return -1
        else:
            return 0
    
    def is_collision(self, A_Vector, B_Vector, C_Vector, P_Vector, V_Vector):
        return ((self.cross_3(A_Vector, B_Vector, P_Vector) == self.cross_3(B_Vector, C_Vector, P_Vector)) and 
                (self.cross_3(B_Vector, C_Vector, P_Vector) == self.cross_3(C_Vector, A_Vector, P_Vector)) and
                ((B_Vector + C_Vector - 2*A_Vector).dot(V_Vector) < 0))
    
    def collision(self, other_Atom):
        P = other_Atom.pos
        V = other_Atom.vel
        if self.is_collision(self.O, self.A, self.B, P, V):
            OT = self.A + self.B - 2*self.O
            other_Atom.vel = other_Atom.vel - 2*OT.dot(other_Atom.vel)*OT/OT.dot(OT)
            
        elif self.is_collision(self.O, self.B, self.C, P, V):
            OT = self.B + self.C - 2*self.O
            other_Atom.vel = other_Atom.vel - 2*OT.dot(other_Atom.vel)*OT/OT.dot(OT)
            
        elif self.is_collision(self.O, self.C, self.D, P, V):
            OT = self.C + self.D - 2*self.O
            other_Atom.vel = other_Atom.vel - 2*OT.dot(other_Atom.vel)*OT/OT.dot(OT)
            
        elif self.is_collision(self.O, self.D, self.A, P, V):
            OT = self.D + self.A - 2*self.O
            other_Atom.vel = other_Atom.vel - 2*OT.dot(other_Atom.vel)*OT/OT.dot(OT)
            
class Element:
    def __init__(self, name, mass, radius, k, c, distance, color):
        self.name = name
        self.mass = mass
        self.radius = radius
        self.color = color
    
    def __str__(self):
        return ('Element(name = ' + name + ', mass = ' + str(self.mass) + 
                ', radius = ' + str(self.radius) + 
                ', color = ' + str(self.color) + ')')
        
class Atom:
    def __init__(self, element, pos, vel = Vector(0, 0)):
        self.element = element
        self.pos = pos
        self.vel = vel
        
    def __str__(self):
        return 'Atom(element = ' + self.element.name + ', pos(' + str(self.pos.x) + ', ' + str(self.pos.y) + '), vel(' + str(self.vel.x) + ', ' + str(self.vel.y) + '))'
    
    def is_collision(self, other_Atom):
        pass
    
    def collision(self, other_Atom):
        d = self.pos - other_Atom.pos
        if not(d.dot(d) == 0):
            m1 = self.element.mass
            v1 = self.vel
            m2 = other_Atom.element.mass
            v2 = other_Atom.vel
            
            v1_ = v1.dot(d)*d/(d.dot(d))
            v2_ = v2.dot(d)*d/(d.dot(d))
            
            if (d.dot(d) < (self.element.radius + other_Atom.element.radius)**2) and (d.dot(v1_-v2_) < 0):
                v1__ = (m1-m2)/(m1+m2)*v1_ + 2*m2/(m1+m2)*v2_
                v2__ = 2*m1/(m1+m2)*v1_ + (m2-m1)/(m1+m2)*v2_
                self.vel = v1 - v1_ + v1__
                other_Atom.vel = v2 - v2_ + v2__
         
class Render:
    def __init__(self, screen, width, height):
        pg.init()
        self.screen = screen
        self.render_vector = Vector(0, height)
        self.render_metric = Tensor(1, 0, 0, -1)
        self.origin_vector = Vector(width/2, height/2)
        
    def rendering_vector(self, vector):
        return self.render_vector + self.render_metric.dot(vector + self.origin_vector)
    
    def text(self, text, font, size, pos, color):
        font_ = pg.font.SysFont(font, size)
        text_ = font_.render(text, True, color)
        render_pos = self.rendering_vector(pos)
        self.screen.blit(text_, (render_pos.x, render_pos.y)) 
    
    def polygon(self, positions, color):
        P_list = []
        for pos in positions:
            P = self.rendering_vector(pos)
            P_list.append([P.x, P.y])
        pg.draw.aalines(self.screen, color, True, P_list, True)
        
    def circle(self, pos, radius, color):
        render_pos = self.rendering_vector(pos)
        pg.draw.circle(self.screen, color, (render_pos.x, render_pos.y), radius)
        
    def wall(self, wall):
        self.polygon([wall.A, wall.B, wall.C, wall.D], wall.color)
        
    def atom(self, atom):
        self.circle(atom.pos, atom.element.radius, atom.element.color)