import pygame as pg
import numpy as np

alive = (255,255,255)
died = (0, 0, 0)

class GameOfLife(object):
    w_width, w_height = 600, 600
    fin = False
    FPS = 60
    screen = pg.display.set_mode((w_width, w_height))

    def __init__(self, matrix, size=25):
        self.size = size
        self.grid = matrix
        self.input_mode = True

     def draw_grid(self):
        cell_size = self.w_width // self.size
        for x in range(self.size):
            for y in range(self.size):
                rect = pg.Rect(x*cell_size, y*cell_size, cell_size, cell_size)
                if self.grid[x][y] == 1:
                    pg.draw.rect(self.screen, alive, rect)
                else:
                    pg.draw.rect(self.screen, died, rect)

    def update_state(self):
        new_grid = np.zeros((self.size, self.size))

        for i in range(self.size):
            for j in range(self.size):
                left = (j-1) % self.size
                right = (j+1) % self.size
                up = (i-1) % self.size
                down = (i+1) % self.size

                neighbors = [self.grid[i, left], self.grid[i, right],
                            self.grid[up, j], self.grid[down, j],
                            self.grid[up, left], self.grid[up, right],
                            self.grid[down, left], self.grid[down, right]]
                total = sum(neighbors)

                if self.grid[i, j] == 1:
                    new_grid[i, j] = 1 if 2 <= total <= 3 else 0
                else:
                    new_grid[i, j] = 1 if total == 3 else 0

        self.grid = new_grid
    
    def run(self):
        pg.init()
        clock = pg.time.Clock()

        while not self.finished:
            clock.tick(self.FPS)
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    self.finished = True

            self.update_state()
            self.draw_grid()

            pg.display.flip()

        pg.quit()

def gen_matrix(volume,size=25):

    grid = np.zeros((size, size))
    indices = np.random.choice(size*size, n, replace=False)
    np.put(grid, indices, 1)

    return grid

volume= 30
matrix = gen_matrix(size)

def main() :
    game = GameOfLife(matrix)
    game.run()

if __name__ == "__main__":
    main()