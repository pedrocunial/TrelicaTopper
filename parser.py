import sys


class Parser:

    def __init__(self, input_file):
        with open(input_file, 'r') as self.fin:
            for i in range(7):
                code = self.fin.readline().strip()
                if code == '*COORDINATES':
                    self.read_coordinates()


    def read_coordinates(self):
        n = self.fin.readline().strip()
        self.coords = []
        for i in range(n):
            line = fin.readline().strip().split()
            self.coords[i][0], coords[i][1] = line[1], line[2]


if __name__ == '__main__':
    input_file = sys.argv[1]
    parser = Parser(input_file)

    print(parser.coords)
