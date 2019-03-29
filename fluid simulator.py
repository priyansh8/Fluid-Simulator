import sys
import pygame
from math import sqrt
from random import random as rand
from math import pi

class Particle(object):
    def __init__(self, pos):
        # Scalars
        self.density = 0

        # Forces
        self.position = pos
        self.velocity = Vec2(0,0)
        self.pressure_force = Vec2(0,0)
        self.viscosity_force = Vec2(0,0)


class ParticleGraphics(object):
    def __init__(self, window_size):
        pygame.init() 
        self.window = pygame.display.set_mode(window_size) 
        self.radius = 5

    def draw(self, particles, H, scale): 
        scale = 500 // scale
        self.window.fill((0,0,0))

        for particle in particles:

            # Color based on pressure just for fun
            if particle.density*50 < 255:
                color = particle.density*50
            else:
                color= 255

            # Area of influence (H)
            pygame.draw.circle(self.window, (color, 0, 255), (
                int(particle.position.x*scale), 
                int(particle.position.y*scale)), 
                H*scale//2, 1)
 
            # Particles
            pygame.draw.circle(self.window, (255, 255, 255), (
                int(particle.position.x*scale), 
                int(particle.position.y*scale)), 
                self.radius)
            
            # Velocity vectors
            pygame.draw.line(self.window, (255, 255, 0), 
                # start
                (int(particle.position.x*scale), 
                 int(particle.position.y*scale)), 
                # end
                (int((particle.position.x+particle.velocity.x)*scale),
                 (int((particle.position.y+particle.velocity.y)*scale))))

        # Enable to step through the sim by entering newlines in the console
        ### getch = input()

        # Draw to screen
        pygame.display.flip() 

graphics = ParticleGraphics((500, 500))

#Num = float, int # Numeric classes

class Vec2(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __add__(self, other):
        return Vec2(self.x+other.x, self.y+other.y)
        
    def __sub__(self, other):
        return Vec2(self.x-other.x, self.y-other.y)

    def __iadd__(self,other):
        self.x += other.x
        self.y += other.y
        return self

    def __eq__(self, other):
        if self.x == other.x and self.y == other.y:
            return True
        else:
            return False

    def __mul__(self, scalar): # vec * scalar
        return Vec2(self.x*scalar, self.y*scalar)

    def __rmul__(self, scalar): # scalar * vec
        return Vec2(self.x*scalar, self.y*scalar)

def length(v):
    return sqrt(v.x**2 + v.y**2)

#MAIN

# Constants
SCALE = 15
MASS = 5 # inversly proportional Particle mass
DENSITY = 1 # Rest density
GRAVITY = Vec2(0, 0.5)
H = 1  # Smoothing cutoff- essentially, particle size
k = 20  # Temperature constant- higher means particle repel more strongly
eta = 1  # Viscosity constant- higher for more viscous


def W(r, h):  # :: Num
    '''
    A weighting function (kernel) for the contribution of each neighbor
    to a particle's density. Forms a nice smooth gradient from the center 
    of a particle to H, where it's 0
    '''
    if 0 < length(r) <= h:
        ret= 315/(64 * pi * h**9) * (h**2 - length(r)**2)**3
    else:
        ret= 0

    # Typecheck
    return ret


def gradient_Wspiky(r, h):  # :: Vec2
    '''
    Gradient ( that is, Vec2(dx, dy) ) of a weighting function for
    a particle's pressure. This weight function is spiky (not flat or
    smooth at x=0) so particles close together repel strongly
    '''
    len_r = length(r)

    if 0 < len_r <= h:
        ret =  -1 * r * (45/(pi * h**6 * len_r)) * (h - len_r)**2
    else:
        ret = Vec2(0, 0)

    return ret
    

def laplacian_W_viscosity(r, h):  # :: Num
    '''
    The laplacian of a weighting function that tends towards infinity when 
    approching 0 (slows down particles moving faster than their neighbors)
    '''
    len_r = length(r)

    if 0 < len_r <= h:
        ret = 45/(2 * pi * h**5) * (1 - len_r/h)
    else:
        ret = 0

    return ret


# Instantiate particles!
##width = 20
##height = 10

particles = []
for x in range(10):
    for y in range(10):
        particles.append(Particle(Vec2(x+1+rand()*0.1, y+5)))


# random distribution
# particles = [Particle(Vec2(rand()*SCALE, rand()*SCALE)) 
#                     for p in range(NUM_PARTICLES)]

time = 0
delta_time = 0.1
while True:

    # Clear everything
    for particle in particles:
        particle.density = DENSITY
        particle.pressure_force = Vec2(0,0)
        particle.viscosity_force = Vec2(0,0)
    
    # Calculate fluid density around each particle
    for particle in particles:
        for neighbor in particles:

            # If particles are close together, density increases
            distance = particle.position - neighbor.position # A vector

            if length(distance) <= H:  # Particles are close enough to matter
                particle.density += MASS * W(distance, H)

    # Calculate forces on each particle based on density
    for particle in particles:
        for neighbor in particles:

            distance = particle.position - neighbor.position
            if length(distance) <= H:
                # Temporary terms used to caclulate forces
                density_p = particle.density
                density_n = neighbor.density
                assert(density_n != 0)  # Dividing by density later

                # Pressure derived 
                pressure_p = k * (density_p - DENSITY)
                pressure_n = k * (density_n - DENSITY)

                # Navier-Stokes equations for pressure and viscosity
                # (ignoring surface tension)
                particle.pressure_force += (-1 *
                        MASS * (pressure_p + pressure_n) / (2 * density_n)
                        * gradient_Wspiky(distance, H))

                particle.viscosity_force += (
                        eta * MASS * (neighbor.velocity - particle.velocity)
                        * (1/density_n) * laplacian_W_viscosity(distance, H))

    # Apply forces to particles- make them move!
    for particle in particles:
        total_force = particle.pressure_force + particle.viscosity_force

        # 'Eulerian' style momentum:

        # Calculate acceleration from forces
        acceleration = total_force * (1/particle.density) \
                * delta_time + GRAVITY

        # Update position and velocity
        particle.velocity += acceleration * delta_time
        particle.position += particle.velocity * delta_time

        # Make sure particles stay in bounds
        # TODO: Better boundary conditions (THESE ARE BAD)
        if particle.position.x >= SCALE - 0.01:
            particle.position.x = SCALE - (0.01 + 0.1*rand())
            particle.velocity.x = 0
        elif particle.position.x < 0.01:
            particle.position.x = 0.01 + 0.1*rand()
            particle.velocity.x = 0

        if particle.position.y >= SCALE - 0.01:
            particle.position.y = SCALE - (0.01+rand()*0.1)
            particle.velocity.y = 0
        elif particle.position.y < 0.01:
            particle.position.y = 0.01 + rand()*0.1
            particle.velocity.y = 0

    graphics.draw(particles, H, SCALE)

    time += delta_time

