# -*- coding: utf-8 -*-
"""
Created on Tue May 17 19:57:28 2022

@author: 39339
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

# Let's make a mini galaxy (9x9)
mini_galaxy = -np.ones((9, 9)) * 2
for i in range(9):
    for j in range(9):
        r2 = (i - 4) * (i - 4) + (j - 4) * (j - 4)
        if r2 >= 1 * 1 and r2 <= 2.5 * 2.5:
            mini_galaxy[i][j] = 0
        elif r2 < 3.5 * 3.5:
            mini_galaxy[i][j] = -1
# sns.heatmap(mini_galaxy)
# plt.show()

# Now let's make a bigger galaxy (100x100)
galaxy = -np.ones((100, 100)) * 2
for i in range(100):
    for j in range(100):
        r2 = (i - 49.5) * (i - 49.5) + (j - 49.5) * (j - 49.5)
        if r2 >= 21 * 21 and r2 <= 36 * 36:
            galaxy[i][j] = 0
        elif r2 < 50 * 50:
            galaxy[i][j] = -1
fig, ax = plt.subplots()
fig.tight_layout()
im = ax.imshow(X=galaxy)
cbar = ax.figure.colorbar(im, ax=ax)
input("Should we continue?")
# For each cell, we want to measure how long it has been since something changed
clock = np.zeros((100, 100))

# Characteristic times - fixed (EVOLUTION)
t01 = 600
t12 = 120
t23 = 20

# Characteristic times - to play with (CATASTROPHIC EVENT)
t10 = 100

t20 = 100
t21 = 100

t30 = 200
t31 = 120
t32 = 40

# Characteristic times - to play with (COLONIZATION)
t03 = 10
t13 = 10
t23_colonization = 10

# Amplitudes - fixed (EVOLUTION)
A01 = 0.003
A12 = 0.003
A23 = 0.003

# Amplitudes - to play with (COLONIZATION)
A03 = 0.5
A13 = 1
A23_colonization = 1

# Amplitudes - to play with (CATASTROPHIC EVENT)
A10 = 0.005

A20 = 0.005
A21 = 0.05

A30 = 0.5
A31 = 0.5
A32 = 0.5


def check_neighbors(galaxy, i, j):
    N = galaxy.shape[0]
    number_of_planets = 0
    number_of_neighbors = 0
    neighbor_planets = []
    for k in [i - 1, i, i + 1]:
        for l in [j - 1, j, j + 1]:
            if k < 0 or k >= N or l < 0 or l >= N:
                continue
            if galaxy[k][l] == -1 or galaxy[k][l] == -2:
                continue
            if k == i and l == j:
                continue
            else:
                number_of_planets += 1
                if galaxy[k][l] == 3:
                    number_of_neighbors += 1
                    neighbor_planets.append((k, l))
    if number_of_neighbors >= 1:
        index = np.random.randint(number_of_neighbors)
        selected_neighbor = neighbor_planets[index]
    else:
        selected_neighbor = None
    return number_of_planets, number_of_neighbors, selected_neighbor


# SIMULATION
galaxy = -np.ones((100, 100)) * 2
for i in range(100):
    for j in range(100):
        r2 = (i - 49.5) * (i - 49.5) + (j - 49.5) * (j - 49.5)
        if r2 >= 21 * 21 and r2 <= 36 * 36:
            galaxy[i][j] = 0
        elif r2 < 50 * 50:
            galaxy[i][j] = -1
world_map = np.copy(galaxy)
civilization_number = 0

for step in range(1000):
    if step % 200 == 0:
        print("Simulation reached", (step // 200), "billion years.")
    # We make a copy of our galaxy's state so that we don't mess it up
    new_galaxy = np.copy(galaxy)
    # We go cell by cell to see what happens there
    for i in range(49 - 36, 50 + 36 + 1):
        for j in range(50 - 36, 50 + 36 + 1):
            # Check if outside the Galactic Habitable Zone
            r2 = (i - 49.5) * (i - 49.5) + (j - 49.5) * (j - 49.5)
            if r2 < 21 * 21 or r2 > 36 * 36:
                # It is outside - we skip this cell
                continue
            planets, neighbors, chosen_neighbor = check_neighbors(new_galaxy, i, j)

            if galaxy[i][j] == 0:
                # NO LIFE ON THIS PLANET
                p1 = A01 * clock[i][j] / t01  # evolution
                p2 = 0
                p3 = (A03 * clock[i][j] / t03) * (neighbors / planets)  # colonization

                p_sum = p1 + p2 + p3
                if p_sum > 1:
                    p1 = p1 / p_sum
                    p2 = p2 / p_sum
                    p3 = p3 / p_sum
                p0 = 1.0 - p_sum

                # Roll the "dice"
                roll = np.random.random()
                if roll < p0:
                    clock[i][j] += 1
                elif roll < p0 + p1:
                    new_galaxy[i][j] = 1
                    clock[i][j] = 0
                elif roll < p0 + p1 + p2:
                    new_galaxy[i][j] = 2
                    clock[i][j] = 0
                else:
                    new_galaxy[i][j] = 3
                    clock[i][j] = 0
                    k, l = chosen_neighbor
                    world_map[i][j] = world_map[k][l]
            elif galaxy[i][j] == 1:
                # SIMPLE LIFE ON THIS PLANET (BACTERIA)
                p0 = A10 * clock[i][j] / t10  # catastrophe
                p2 = A12 * clock[i][j] / t12  # evolution
                p3 = (A13 * clock[i][j] / t13) * (neighbors / planets)  # colonization

                p_sum = p0 + p2 + p3
                if p_sum > 1:
                    p0 = p0 / p_sum
                    p2 = p2 / p_sum
                    p3 = p3 / p_sum
                p1 = 1.0 - p_sum

                # Roll the "dice"
                roll = np.random.random()
                if roll < p0:
                    new_galaxy[i][j] = 0
                    clock[i][j] = 0
                elif roll < p0 + p1:
                    clock[i][j] += 1
                elif roll < p0 + p1 + p2:
                    new_galaxy[i][j] = 2
                    clock[i][j] = 0
                else:
                    new_galaxy[i][j] = 3
                    clock[i][j] = 0
                    k, l = chosen_neighbor
                    world_map[i][j] = world_map[k][l]
            elif galaxy[i][j] == 2:
                # COMPLEX LIFE ON THIS PLANET (PLANTS AND ANIMALS)
                p0 = A20 * clock[i][j] / t20  # catastrophe
                p1 = A21 * clock[i][j] / t21  # evolution
                p3_colonization = (
                    A23_colonization * clock[i][j] / t23_colonization
                ) * (
                    neighbors / planets
                )  # colonization
                p3_evolution = A23 * clock[i][j] / t23  # evolution

                p_sum = p0 + p1 + p3_colonization + p3_evolution
                if p_sum > 1:
                    p0 = p0 / p_sum
                    p1 = p1 / p_sum
                    p3_colonization = p3_colonization / p_sum
                    p3_evolution = p3_evolution / p_sum
                p2 = 1.0 - p_sum

                # Roll the "dice"
                roll = np.random.random()
                if roll < p0:
                    new_galaxy[i][j] = 0
                    clock[i][j] = 0
                elif roll < p0 + p1:
                    new_galaxy[i][j] = 1
                    clock[i][j] = 0
                elif roll < p0 + p1 + p2:
                    clock[i][j] += 1
                elif roll < p0 + p1 + p2 + p3_colonization:
                    new_galaxy[i][j] = 3
                    clock[i][j] = 0
                    k, l = chosen_neighbor
                    world_map[i][j] = world_map[k][l]
                else:
                    new_galaxy[i][j] = 3
                    clock[i][j] = 0
                    civilization_number += 1
                    world_map[i][j] = civilization_number
            elif new_galaxy[i][j] == 3:
                # TECHNOLOGICALLY ADVANCED CIVILIZATION
                p0 = A30 * clock[i][j] / t30  # catastrophe
                p1 = A31 * clock[i][j] / t31  # catastrophe
                p2 = A32 * clock[i][j] / t32  # catastrophe

                p_sum = p0 + p1 + p2
                if p_sum > 1:
                    p0 = p0 / p_sum
                    p1 = p1 / p_sum
                    p2 = p2 / p_sum
                p3 = 1.0 - p_sum

                # Roll the "dice"
                roll = np.random.random()
                if roll < p0:
                    new_galaxy[i][j] = 0
                    clock[i][j] = 0
                    world_map[i][j] = 0
                elif roll < p0 + p1:
                    new_galaxy[i][j] = 1
                    clock[i][j] = 0
                    world_map[i][j] = 0
                elif roll < p0 + p1 + p2:
                    new_galaxy[i][j] = 2
                    clock[i][j] = 0
                    world_map[i][j] = 0
                else:
                    clock[i][j] += 1
    galaxy = np.copy(new_galaxy)
print("Simulation finished at 5 billion years")

plt.clf()
plt.cla()
plt.figure(figsize=(10, 8))
fig, ax = plt.subplots()
fig.tight_layout()
im = ax.imshow(X=galaxy)
cbar = ax.figure.colorbar(im, ax=ax)
plt.show()

input("Do you want to see the world map?")

plt.clf()
plt.cla()
plt.figure(figsize=(10, 8))
fig, ax = plt.subplots()
fig.tight_layout()
im = ax.imshow(X=galaxy)
cbar = ax.figure.colorbar(im, ax=ax)
plt.show()
