# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 14:16:50 2021

@author: junob
"""
import math as m
from vectortools import *
import pygame as pg

def cross_3(A_Vector, B_Vector, C_Vector):
    AB = A_Vector - B_Vector
    AC = A_Vector - C_Vector
    if AB.cross(AC) > 0:
        return 1
    elif AB.cross(AC) < 0:
        return -1
    else:
        return 0

def is_in_triangle(A_Vector, B_Vector, C_Vector, P_Vector):
    return ((cross_3(A_Vector, B_Vector, P_Vector) == cross_3(B_Vector, C_Vector, P_Vector)) and 
            (cross_3(B_Vector, C_Vector, P_Vector) == cross_3(C_Vector, A_Vector, P_Vector)))
    
class Wall:
    def __init__(self, width, height, theta, pos, color):
        self.width = width
        self.height = height
        self.theta = theta
        self.pos = pos
        self.color = color
        self.O = pos + SO2(theta).dot(Vector(width/2, height/2))
        
        A = pos
        B = pos + SO2(theta).dot(Vector(width, 0))
        C = pos + SO2(theta).dot(Vector(width, height))
        D = pos + SO2(theta).dot(Vector(0, height))
        self.P = [A, B, C, D]
        
    def __str__(self):
        return ('Wall(width = ' + str(self.width) + ', height = ' + str(self.height) + 
                ', theta = ' + str(self.theta) + 
                ', pos(' + str(self.pos.x) + ', ' + str(self.pos.y) + ')' + 
                ', color = ' + str(self.color) + ')')
    
    def is_collision(self, other_Atom):
        P = other_Atom.pos
        V = other_Atom.vel
        for i in range(4):
            if is_in_triangle(self.O, self.P[i-1], self.P[i], P):
                OT = self.P[i-1] + self.P[i] - 2*self.O
                if (OT.dot(V) < 0):
                    return True
        return False
            
class Element:
    def __init__(self, name, mass, radius, color):
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
    
    def is_collision(self, other):
        if isinstance(other, Atom):
            d = self.pos - other.pos
            if not(d.dot(d) == 0):
                v1 = self.vel
                v2 = other.vel
                
                v1_ = v1.dot(d)*d/(d.dot(d))
                v2_ = v2.dot(d)*d/(d.dot(d))
                
                return (d.dot(d) < (self.element.radius + other.element.radius)**2) and (d.dot(v1_-v2_) < 0)
            else:
                return False
        
        elif isinstance(other, Wall):
            P = self.pos
            V = self.vel
            for i in range(4):
                if is_in_triangle(other.O, other.P[i-1], other.P[i], P):
                    OT = other.P[i-1] + other.P[i] - 2*other.O
                    if (OT.dot(V) < 0):
                        return True
            return False
                        
    def collision(self, other):
        if isinstance(other, Atom):
            d = self.pos - other.pos
            if not(d.dot(d) == 0):
                m1 = self.element.mass
                v1 = self.vel
                m2 = other.element.mass
                v2 = other.vel
                
                v1_ = v1.dot(d)*d/(d.dot(d))
                v2_ = v2.dot(d)*d/(d.dot(d))
                
                if (d.dot(d) < (self.element.radius + other.element.radius)**2) and (d.dot(v1_-v2_) < 0):
                    v1__ = (m1-m2)/(m1+m2)*v1_ + 2*m2/(m1+m2)*v2_
                    v2__ = 2*m1/(m1+m2)*v1_ + (m2-m1)/(m1+m2)*v2_
                    self.vel = v1 - v1_ + v1__
                    other.vel = v2 - v2_ + v2__
            
        elif isinstance(other, Wall):
            P = self.pos
            V = self.vel
            for i in range(4):
                if is_in_triangle(other.O, other.P[i-1], other.P[i], P):
                    OT = other.P[i-1] + other.P[i] - 2*other.O
                    if (OT.dot(V) < 0):
                        self.vel = V - 2*OT.dot(V)*OT/OT.dot(OT)
                        
class World:
    def __init__(self, t, atoms, walls, gravity):
        self.t = t
        self.atoms = atoms
        self.walls = walls
        self.gravity = gravity
        
    def __str__(self):
        return ('World(t = ' + self.t + ', atoms = ' + str(self.atoms) + 
                ', walls = ' + str(self.walls) + ', gravity = ' + str(self.gravity) + ')')
    
class Render:
    def __init__(self, screen, width, height):
        pg.init()
        self.screen = screen
        self.width = width
        self.height = height
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
        self.polygon(wall.P, wall.color)
        
    def atom(self, atom):
        self.circle(atom.pos, atom.element.radius, atom.element.color)
        
class Simulator:
    def __init__(self, dt, world, render):
        self.dt = dt
        self.world = world
        self.render = render
        self.count = 0
        
    def clock(self):
        self.world.t = self.world.t + self.dt
        return self.world.t
    
    def draw_background(self, color):
        self.render.screen.fill(color)
    
    def draw_grid(self, unit_size):
        for x in range(0, int(self.render.width/2), unit_size):
            grey = (200, 200, 200)
            pg.draw.line(self.render.screen, grey, (x + self.render.width/2, 0), (x + self.render.width/2, self.render.height))
            pg.draw.line(self.render.screen, grey, (-x + self.render.width/2, 0), (-x + self.render.width/2, self.render.height))
        for y in range(0, int(self.render.height/2), unit_size):
            pg.draw.line(self.render.screen, grey, (0, y + self.render.height/2), (self.render.width, y + self.render.height/2))
            pg.draw.line(self.render.screen, grey, (0, -y + self.render.height/2), (self.render.width, -y + self.render.height/2))
        
    def draw_wall(self):
        for wall in self.world.walls:
            self.render.wall(wall)
    
    def draw_atom(self):
        for atom in self.world.atoms:
            self.render.atom(atom)
            
    def atom_atom_collision(self):
        for atom in self.world.atoms:
            for other_atom in self.world.atoms:
                atom.collision(other_atom)
                
    def atom_wall_collision(self):
        for atom in self.world.atoms:
            for wall in self.world.walls:
                atom.collision(wall)
                
    def main(self):
        x_ = []
        v_ = []
        for atom in self.world.atoms:
            new_v = atom.vel + self.world.gravity*self.dt
            v_.append(new_v)
            x_.append(atom.pos + new_v*self.dt)
            
        count = 0
        for atom in self.world.atoms:
            atom.pos = x_[count]
            atom.vel = v_[count]
            count = count + 1
 
    def save_screen(self, directory):
        self.count += 1
        img = directory + '/%08d.png' % (self.count)
        pg.image.save(self.render.screen, img)
    
if __name__ == '__main__':
    width = 1000
    height = 800
    
    screen = pg.display.set_mode((width, height))
    render = Render(screen, width, height)
    clock = pg.time.Clock() 
    
    black = pg.Color('black')
    white = pg.Color('white')
    red = pg.Color('red')
    green = pg.Color('green')
    blue = pg.Color('blue')
    
    wall1 = Wall(1000, 50, 0, Vector(-500, -400), blue)
    wall2 = Wall(50, 800, 0, Vector(-500, -400), blue)
    wall3 = Wall(50, 800, 0, Vector(450,-400), blue)
    wall4 = Wall(1000, 50, 0, Vector(-500, 350), blue)
    wall5 = Wall(100, 50, m.pi/4, Vector(-300, 0), blue)
    
    e1 = Element(name = 'Heilum', mass = 1, radius = 10, color = red)
    atom1 = Atom(e1, Vector(-200, 0), Vector(50, 0))
    atom2 = Atom(e1, Vector(0, 0))
    atom3 = Atom(e1, Vector(25, -10))
    atom4 = Atom(e1, Vector(25, 10))
    atom5 = Atom(e1, Vector(50, -20))
    atom6 = Atom(e1, Vector(50, 0))
    atom7 = Atom(e1, Vector(50, 20))
    
    walls = [wall1, wall2, wall3, wall4, wall5]
    atoms = [atom1, atom2, atom3, atom4, atom5, atom6, atom7]
    
    gravity = Vector(0, 0)
    
    world = World(0, atoms, walls, gravity)
    
    simulator = Simulator(0.01, world, render)
    
    while True:
        t = simulator.clock()
        simulator.draw_background(white)
        simulator.draw_grid(100)
        simulator.draw_wall()
        simulator.atom_wall_collision()
        simulator.atom_atom_collision()
        simulator.main()
        simulator.draw_atom()
        
        render.text('pos = (%.2f, %.2f)'%(atom1.pos.x, atom1.pos.y) , None, 30, Vector(atom1.pos.x -100, atom1.pos.y - 30), black)
        render.text('vel = (%.2f, %.2f)'%(atom1.vel.x, atom1.vel.y) , None, 30, Vector(atom1.pos.x -100, atom1.pos.y - 50), black)
        
        render.text('pos = (%.2f, %.2f)'%(atom7.pos.x, atom7.pos.y) , None, 30, Vector(atom7.pos.x -100, atom7.pos.y - 30), blue)
        render.text('vel = (%.2f, %.2f)'%(atom7.vel.x, atom7.vel.y) , None, 30, Vector(atom7.pos.x -100, atom7.pos.y - 50), blue)
        
        for event in pg.event.get():
            if event.type == pg.QUIT:
                sys.exit()
        clock.tick(100)
        pg.display.update()
        
        simulator.save_screen('images/pocket_ball_demo')
