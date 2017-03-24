import sys


class Parser:

    INCI_ARG_COUNT = 3


    def __init__(self, input_file):
        with open(input_file, 'r') as self.fin:
            for line in self.fin:
                line = line.strip()
                if line == '*COORDINATES':
                    self.read_coordinates()
                elif line == '*INCIDENCES':
                    self.read_incidences()
                elif line == '*MATERIALS':
                    self.read_materials()
                elif line == '*GEOMETRIC_PROPERTIES':
                    self.read_properties()
                elif line == '*BCNODES':
                    self.read_bcnodes()
                elif line == '*LOADS':
                    self.read_loads()



    def read_coordinates(self):
        n = int(self.fin.readline().strip())
        self.coords = []
        for i in range(n):
            line = self.fin.readline().strip().split()
            tmp = [0, 0]
            tmp[0], tmp[1] = float(line[1]), float(line[2])
            self.coords.append(tmp)


    def read_incidences(self):
        self.inci = []
        line = self.fin.readline().strip().split()
        while len(line) == self.INCI_ARG_COUNT:
            tmp = [0, 0]
            tmp[0], tmp[1] = int(line[1]), int(line[2])
            self.inci.append(tmp)
            line = self.fin.readline().strip().split()


    def read_materials(self):
        self.materials = []
        n = int(self.fin.readline().strip())
        for i in range(n):
            line = self.fin.readline().strip().split()
            self.materials.append([float(line[0])])


    def read_properties(self):
        self.prop = []
        n = int(self.fin.readline().strip())
        for i in range(n):
            line = self.fin.readline().strip().split()
            self.prop.append([float(line[0])])


    def read_bcnodes(self):
        self.bcnodes = []
        n = int(self.fin.readline().strip())
        for i in range(n):
            line = self.fin.readline().strip().split()
            self.bcnodes.append([int(line[0]), int(line[1])])


    def read_loads(self):
        self.loads = []
        n = int(self.fin.readline().strip())

        # this can be optimized
        for i in range(len(self.coords) * 2):
            self.loads.append(0)

        for i in range(n):
            line = self.fin.readline().strip().split()
            self.loads[(int(line[0]) - 1) * 2 + int(line[1]) - 1] = int(line[2])



if __name__ == '__main__':
    input_file = sys.argv[1]
    parser = Parser(input_file)

    print(parser.coords)
    print(parser.inci)
    print(parser.materials)
    print(parser.prop)
    print(parser.bcnodes)
    print(parser.loads)
