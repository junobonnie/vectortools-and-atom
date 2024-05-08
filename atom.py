# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 14:16:50 2021

@author: junob
"""
import math as m
from vectortools import *
import pygame as pg
import sys
import h5py

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
        return ('Element(name = ' + self.name + ', mass = ' + str(self.mass) +
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
            if not self == other:
                d = self.pos - other.pos
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

    def collision(self, other, dt):
        if isinstance(other, Atom):
            if not self == other:
                d = self.pos - other.pos
                m1 = self.element.mass
                v1 = self.vel
                m2 = other.element.mass
                v2 = other.vel

                v1_ = v1.dot(d)*d/(d.dot(d))
                v2_ = v2.dot(d)*d/(d.dot(d))

                if (d.dot(d) < (self.element.radius + other.element.radius)**2) and (d.dot(v1_-v2_) < 0):
                    self.pos -= self.vel*dt
                    other.pos -= other.vel*dt
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
                        self.pos -= self.vel*dt
                        self.vel = V - 2*OT.dot(V)*OT/OT.dot(OT)

    def fusion(self, other_Atom):
        new_Atom = None
        if not self == other_Atom:
            d = self.pos - other_Atom.pos
            if (d.dot(d) < (self.element.radius + other_Atom.element.radius)**2):
                new_element = Element(name = 'New atom', mass = self.element.mass + other_Atom.element.mass, 
                                      radius = m.sqrt(self.element.radius**2 + other_Atom.element.radius**2),
                                      color = self.element.color + other_Atom.element.color)
                new_Atom = Atom(element = new_element, 
                                pos = (self.element.mass*self.pos + other_Atom.element.mass*other_Atom.pos)/(self.element.mass + other_Atom.element.mass),
                                vel = (self.element.mass*self.vel + other_Atom.element.mass*other_Atom.vel)/(self.element.mass + other_Atom.element.mass))
        return new_Atom
    
class World:
    def __init__(self, t, atoms, walls, gravity):
        self.t = t
        self.atoms = atoms
        self.walls = walls
        self.gravity = gravity

    def __str__(self):
        return ('World(t = ' + str(self.t) + ', atoms = ' + str(self.atoms) +
                ', walls = ' + str(self.walls) + ', gravity = ' + str(self.gravity) + ')')

class Render:
    def __init__(self, screen, width, height):
        if not screen == None:
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
    def __init__(self, dt, world, render, grid_size = 100):
        self.dt = dt
        self.world = world
        self.render = render
        self.count_screen = 0
        self.count_snapshot = 0
        self.grid_size = grid_size
        self.grid = None

    def clock(self):
        self.world.t = self.world.t + self.dt
        return self.world.t

    def draw_background(self, color):
        self.render.screen.fill(color)

    def draw_grid(self, unit_size):
        grey = (200, 200, 200)
        for x in range(0, int(self.render.width/2), unit_size):
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
            if ((-self.render.width/2 < atom.pos.x < self.render.width/2) and 
                    (-self.render.height/2 < atom.pos.y < self.render.height/2)):
                self.render.atom(atom)

    def make_grid(self):
        nx = int(self.render.width//self.grid_size+1)
        ny = int(self.render.height//self.grid_size+1)
        grid = [[] for i in range(nx*ny)]
        for atom in self.world.atoms:
            i = int((self.render.width/2 + atom.pos.x)//self.grid_size)
            j = int((self.render.height/2 + atom.pos.y)//self.grid_size)
            if (0 <= i < nx) and (0 <= j < ny):
                grid[i+nx*j].append(atom)
        self.grid = grid
    
    def get_near_atoms(self, atom):
        nx = int(self.render.width//self.grid_size+1)
        ny = int(self.render.height//self.grid_size+1)
        i = int((self.render.width/2 + atom.pos.x)//self.grid_size)
        j = int((self.render.height/2 + atom.pos.y)//self.grid_size)
        atoms = []
        for i_ in (i-1, i, i+1):
            for j_ in (j-1, j, j+1):
                if (0 <= i_ < nx) and (0 <= j_ < ny):
                    atoms += self.grid[i_+nx*j_]
        return atoms

    def atom_atom_collision(self):
        self.make_grid()
        for atom in self.world.atoms:
            atoms = self.get_near_atoms(atom)
            for other_atom in atoms:
                #self.render.polygon([atom.pos, other_atom.pos], red)
                atom.collision(other_atom, self.dt)

    def atom_wall_collision(self):
        for atom in self.world.atoms:
            for wall in self.world.walls:
                atom.collision(wall, self.dt)
    
    def atom_atom_fusion(self):
        while True:
            for atom in self.world.atoms[:]:
                for other_atom in self.world.atoms[:]:
                    new_atom = atom.fusion(other_atom)
                    if not new_atom == None:
                        self.world.atoms.remove(atom)
                        self.world.atoms.remove(other_atom)
                        self.world.atoms.append(new_atom)
                        break
                if not new_atom == None:
                    break
            if new_atom == None:
                break
                
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

    def save_screen(self, directory, skip_number = 0):
        if self.count_screen%(skip_number+1) == 0:
            img = directory + '/%08d.png' % (self.count_screen)
            pg.image.save(self.render.screen, img)
        self.count_screen += 1
        
    # def save_snapshot(self, directory, skip_number = 0):
    #     if self.count_snapshot%(skip_number+1) == 0:
    #         snapshot = directory + '/snapshot_%08d.txt' % (self.count_snapshot)
    #         with open(snapshot, "w") as f:
    #             walls_info = ''
    #             count = 0
    #             for wall in self.world.walls:
    #                 count += 1
    #                 walls_info = walls_info + 'wall' + str(count) + '{ width:' + str(wall.width) + ', height:' + str(wall.height) + ', theta:' + str(wall.theta) + ', pos:' + str(wall.pos) + ', color:' + str(wall.color) + ' }, '
                    
    #             atoms_info = ''
    #             count = 0
    #             for atom in self.world.atoms:
    #                 count += 1
    #                 atoms_info = atoms_info + 'atom' + str(count) + '{ element{ name:' + atom.element.name + ', mass:' + str(atom.element.mass) + ', radius:' + str(atom.element.radius) + ', color:' + str(atom.element.color) + ' }' + ', pos:' + str(atom.pos) + ', vel:' + str(atom.vel) + ' }, '
                    
    #             f.write('world{ t:' + str(self.world.t) + ', gravity:' + str(self.world.gravity) + ', walls{ ' + walls_info + ' }' + ', atoms{ ' + atoms_info + ' }' + ' }')
    #     self.count_snapshot += 1
        
    def save_snapshot(self, directory, skip_number = 0):
        if self.count_snapshot%(skip_number+1) == 0:
            snapshot = directory + '/snapshot_%08d.hdf5' % (self.count_snapshot)
            with h5py.File(snapshot, 'w') as f:
                f.attrs['count_snapshot'] = self.count_snapshot
                world = f.create_group('world')
                world.attrs['t'] = self.world.t
                world.attrs['gravity'] = self.world.gravity.list()
                atoms = world.create_group('atoms')
                N = len(self.world.atoms)
                element = [0]*N
                mass = [0]*N
                radius = [0]*N
                color = [0]*N
                pos = [0]*N
                vel = [0]*N
                count = 0
                for atom in self.world.atoms:
                    element[count] = atom.element.name
                    mass[count] = atom.element.mass
                    radius[count] = atom.element.radius
                    color[count] = (atom.element.color.r, atom.element.color.g, atom.element.color.b, atom.element.color.a)
                    pos[count] = atom.pos.list()
                    vel[count] = atom.vel.list()
                    count += 1
                atoms.create_dataset('element', data = element)
                atoms.create_dataset('mass', data = mass)
                atoms.create_dataset('radius', data = radius)
                atoms.create_dataset('color', data = color)
                atoms.create_dataset('pos', data = pos)
                atoms.create_dataset('vel', data = vel)
                walls = world.create_group('walls')
                N = len(self.world.walls)
                width = [0]*N
                height = [0]*N
                theta = [0]*N
                pos = [0]*N
                color = [0]*N
                count = 0
                for wall in self.world.walls:
                    width[count] = wall.width
                    height[count] = wall.height
                    theta[count] = wall.theta
                    pos[count] = wall.pos.list()
                    color[count] = (wall.color.r, wall.color.g, wall.color.b, wall.color.a)
                    count += 1
                walls.create_dataset('width', data = width)
                walls.create_dataset('height', data = height)
                walls.create_dataset('theta', data = theta)
                walls.create_dataset('pos', data = pos)
                walls.create_dataset('color', data = color)
        self.count_snapshot += 1

    def list_to_vector(self, list):
        return Vector(float(list[0]), float(list[1]))

    def load_snapshot(self, snapshot_file):
        with h5py.File(snapshot_file, 'r') as f:
            self.count_snapshot = int(f.attrs['count_snapshot'])
            world = f['world']
            t = float(world.attrs['t'])
            gravity = self.list_to_vector(world.attrs['gravity'])
            element_ = world['atoms']['element']
            mass_ = world['atoms']['mass']
            radius_ = world['atoms']['radius']
            color_ = world['atoms']['color']
            pos_ = world['atoms']['pos']
            vel_ = world['atoms']['vel']
            N = len(element_)
            atoms = [0]*N
            for i in range(N):
                element = Element(element_[i], float(mass_[i]), float(radius_[i]), pg.Color(color_[i]))
                pos = self.list_to_vector(pos_[i])
                vel = self.list_to_vector(vel_[i])
                atoms[i] = Atom(element, pos, vel)
            width_ = world['walls']['width']
            height_ = world['walls']['height']
            theta_ = world['walls']['theta']
            pos_ = world['walls']['pos']
            color_ = world['walls']['color']
            N = len(width_)
            walls = [0]*N
            for i in range(N):
                walls[i] = Wall(float(width_[i]), float(height_[i]), float(theta_[i]), self.list_to_vector(pos_[i]), pg.Color(color_[i]))
            self.world = World(t, atoms, walls, gravity)
        
    # def load_snapshot(self, snapshot_file):
    #     with open(snapshot_file, "r") as f:
    #         snapshot = f.read()
    #         snapshot = snapshot.replace('world{ t:', '').replace(', gravity:', '#').replace(', walls', '#').replace('}, atoms', '#').replace(' },  } }', '') 
    #         snapshot = snapshot.split('#')
    #         t = float(snapshot[0])
    #         gravity = eval(snapshot[1])
    #         walls_raw = snapshot[2]
    #         walls_raw = walls_raw.replace('{ wall', '').split('wall')
    #         walls = []
    #         for wall in walls_raw:
    #             wall = wall.replace('width:', '#').replace(', height:', '#').replace(', theta:', '#').replace(', pos:', '#').replace(', color:', '#').replace(' }, ', '').replace(' }', '')
    #             wall = wall.split('#')
    #             try:
    #                 walls.append(Wall(float(wall[1]), float(wall[2]), float(wall[3]), eval(wall[4]), eval(wall[5])))
    #             except:
    #                 pass
                
    #         atoms_raw = snapshot[3]
    #         atoms_raw = atoms_raw.replace('{ atom', '').split('atom')
    #         atoms = []
    #         for atom in atoms_raw:
    #             atom = atom.replace('element{ ', '#').replace(' }, pos:', '#').replace(', vel:', '#').replace(' }, ', '')
    #             atom = atom.split('#')
    #             try:
    #                 element_raw = atom[1]
    #                 element_raw = element_raw.replace('name:', '#').replace(', mass:', '#').replace(', radius:', '#').replace(', color:', '#')
    #                 element_raw = element_raw.split('#')
    #                 atoms.append(Atom(Element(element_raw[1], float(element_raw[2]), float(element_raw[3]), eval(element_raw[4])), eval(atom[2]), eval(atom[3])))
    #             except:
    #                 pass
    #     self.world = World(t, atoms, walls, gravity)
                
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

    e1 = Element(name = 'Helium', mass = 1, radius = 10, color = red)
    atom1 = Atom(e1, Vector(-200, 0), Vector(50, 0))
    atom2 = Atom(e1, Vector(0, 0))
    atom3 = Atom(e1, Vector(25, -10))
    atom4 = Atom(e1, Vector(25, 10))
    atom5 = Atom(e1, Vector(50, -20))
    atom6 = Atom(e1, Vector(50, 0))
    atom7 = Atom(e1, Vector(50, 20))
    
    walls = [wall1, wall2, wall3, wall4, wall5]
    atoms = [atom1, atom2, atom3, atom4, atom5, atom6, atom7]

    gravity = Vector(0, -10)*0

    world = World(0, atoms, walls, gravity)

    simulator = Simulator(0.01, world, render)
    simulator.load_snapshot('snapshots/pocket_ball_demo/snapshot_00000300.hdf5')

    atom1 = simulator.world.atoms[0]
    atom7 = simulator.world.atoms[6]
    while True:
        t = simulator.clock()
        simulator.draw_background(white)
        simulator.draw_grid(100)
        simulator.draw_wall()
        simulator.atom_wall_collision()
        simulator.atom_atom_collision()
        #simulator.atom_atom_fusion()
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
        
        #simulator.save_screen('images/pocket_ball_demo')
        simulator.save_snapshot('snapshots/pocket_ball_demo', 99)
