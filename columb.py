import math
import matplotlib.pyplot as plt
import os

cmd = 'mode 110,50'
os.system(cmd)

clear = lambda: os.system('cls')

#Global Variables
DATA = [50, 0.5, 75, 12, 10, 20, 10, 5, 30, 15, 15, 3, [(10,1), (20,4)]]
INPUT_MSGS = [("int", "\nselect the number of wedges: "),
              ("float", "\nselect the desired inclined step width (m): " ),
              ("float", "\nselect the desired wall angle (deg): " ),
              ("float", "\nselect the desired backfill angle (deg): " ),
              ("float", "\nselect the desired vertical wall height (m): " ),
              ("float", "\nselect the desired soil saturated density (kN/m3): " ),
              ("float", "\nselect the desired cohesion (shear param) (kPa): " ),
              ("float", "\nselect the desired adhesion (interface param) (kPa): " ),
              ("float", "\nselect the desired internal friction angle (shear param) (deg): " ),
              ("float", "\nselect the desired wall friction angle (interface param) (deg): " ),
              ("float", "\nselect the uniform surcharge value (0 if no surcharge) (kPa): " ),
              ("float", "\nselect the desired GWT level from base (-1 if no GWT) (m): " ),
              ("line", "\nselect consecutively the desired line load(s) and inclined distance(s)\n"
                        "from R.Wall top separated by a spaces (0 0 for no line load) (kN/m'): ")
              ]


#Classes
class FailureWedge:
    def __init__(self, density, coordinates, GWT_level=-1, water_density=10):
        self.coordinates = list(coordinates)
        self.density = density
        self.water_density = water_density
        self.sat_coordinates = []
        self.unsat_coordinates = []

        self.lines = self.gen_lines()
        self.resolve_GWT(GWT_level)
        self.weight = self.calculate_weight()

    def is_homogenous(self, GWT_level):
        if GWT_level == -1: return True
        return False
    def gen_lines(self):
        lines = []  # list of tuples. each tuple has coordinate indicies in self.coordinates for a line
        l = len(self.coordinates)
        for iter in range(l - 1):
            lines.append((iter, iter + 1))
        lines.append((l - 1, 0))

        return lines
    def set_sat_coordinates(self, cor):
        self.sat_coordinates = cor
    def set_unsat_coordinates(self, cor):
        self.unsat_coordinates = cor
    def resolve_GWT(self, GWT_level):
        if self.is_homogenous(GWT_level):
            self.set_unsat_coordinates(self.coordinates)
            self.set_sat_coordinates([])
            return


        intersection_points_with_GWT = []
        line2 = [(0, GWT_level), (1, GWT_level)]

        for line in self.lines:
            line1 = [self.coordinates[line[0]], self.coordinates[line[1]]]
            point = intersection_point(line1,line2)
            if point not in intersection_points_with_GWT: intersection_points_with_GWT.append(point)
            else: intersection_points_with_GWT.append(None)

        indices_of_intersection = [i for i, x in enumerate(intersection_points_with_GWT) if x != None]
        indices_of_intersection.reverse()


        # inserting intersection points in self.coordinates
        for indi in indices_of_intersection:
           self.coordinates.insert(indi+1, intersection_points_with_GWT[indi])

        indices_of_intersection.reverse()
        #updated indicies for GWT intersection points in self.coordinates
        for i in range(len(indices_of_intersection)):
            indices_of_intersection[i] += (i+1)

        #update sat and unsat coordinates

        self.set_unsat_coordinates(self.coordinates[:indices_of_intersection[0]+1] + self.coordinates[indices_of_intersection[1]:])
        self.set_sat_coordinates(self.coordinates[indices_of_intersection[0]:indices_of_intersection[1]+1])


    def get_wedge_cg(self, uniform_surcharge, dista, init_soil_coor, backfill_angle):

        surch_coor = (init_soil_coor[0]+ dista*math.cos(math.radians(backfill_angle))/2,
                      init_soil_coor[1]+ dista*math.sin(math.radians(backfill_angle))/2)
        sur_line = uniform_surcharge * dista

        sat_cg = get_cg(self.sat_coordinates)
        unsat_cg = get_cg(self.unsat_coordinates)
        sub_density = self.density-self.water_density

        if len(self.sat_coordinates) == 0: sub_density = 0

        return ((surch_coor[0]*sur_line + unsat_cg[0]*self.density + sat_cg[0]*(sub_density))/(sur_line+self.density+sub_density),
                (surch_coor[1]*sur_line + unsat_cg[1]*self.density + sat_cg[1]*(sub_density))/(sur_line+self.density+sub_density))

    def calculate_weight(self):
        return area_calculation(self.sat_coordinates) * (self.density - self.water_density) + area_calculation(self.unsat_coordinates) * self.density

###################################################

#Helper Functions
def area_calculation(list_of_coordinates):
    """takes a list of 2D-coordinates(tuples) and reutrns area"""
    if len(list_of_coordinates) == 0: return 0
    area = 0
    l = len(list_of_coordinates)
    for iter in range(l-1):
        area = area + list_of_coordinates[iter][0]*list_of_coordinates[iter+1][1] - list_of_coordinates[iter][1]*list_of_coordinates[iter+1][0]
    area = area + list_of_coordinates[l-1][0]*list_of_coordinates[0][1] - list_of_coordinates[l-1][1]*list_of_coordinates[0][0]
    return abs(area/2)
def get_cg(list_of_coordinates):
    """gets CG of list of coordinates"""
    l = len(list_of_coordinates)
    if l == 0 : return [0, 0]
    s = [0.0, 0.0]
    for coor in list_of_coordinates:
        s[0] += coor[0]
        s[1] += coor[1]

    return(s[0]/l, s[1]/l)
def distance(point1, point2):
    """returns distance between two points"""
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)
def isBetween(middle_point, line):
    """checks if a point lies on a st. line between two other points"""
    if distance(line[0], middle_point) <= distance(line[0], line[1]): return True
    return False
def intersection_point(line1, line2):
    """takes lines as a list of two 2D coordinates(tuples) and returns the coordinate of their point of intersection.
    if no intersection is found returns None"""

    # slopes

    if (line1[0][0] - line1[1][0]) != 0:
        m1 = (line1[0][1] - line1[1][1])/(line1[0][0] - line1[1][0])
        param1 = [m1, -1, line1[0][1] - line1[0][0] * m1]
    else:
        m1 = float('inf')
        param1 = [1, 0, -1 * line1[0][0]]

    if (line2[0][0] - line2[1][0]) != 0:
        m2 = (line2[0][1] - line2[1][1])/(line2[0][0] - line2[1][0])
        param2 = [m2, -1, line2[0][1] - line2[0][0]*m2]
    else:
        m2 = float('inf')
        param2 = [1, 0, -1 * line2[0][0]]

    # return None if parallel
    if m1 == m2: return None



    solution =  ((param1[1]*param2[2] - param2[1]*param1[2])/(param1[0]*param2[1] - param2[0]*param1[1]),
            (param1[2]*param2[0] - param2[2]*param1[0])/(param1[0]*param2[1] - param2[0]*param1[1]))

    # checks if intersection point lies inside line segment
    if not isBetween(solution, line1): return None

    return solution
def line_angle(point1, point2):
    """calculates the slope of a line"""
    if (point1[0]-point2[0]) == 0: return float('inf')
    slope = (point2[1]-point1[1]) / (point2[0]-point1[0])
    return math.degrees(math.atan(slope))
def format_lines(inputted):
    input_line_loads = [float(x) for x in inputted.split()]
    line_loads = []
    for i in range(0, len(input_line_loads), 2):
        line_loads.append(tuple([input_line_loads[i], input_line_loads[i + 1]]))
    return line_loads
def get_line_intersection_with_backfill(wall_coordinates, angle, point):
    arb_point = (point[0] - 1000 * math.cos(math.radians(angle)), point[1] - 1000 * math.sin(math.radians(angle)))
    return intersection_point(wall_coordinates, [point,arb_point])

def get_p_active_point(wall_coordinates, angle, cg):
    arb_point = (cg[0] - 1000 * math.cos(math.radians(angle)), cg[1] - 1000 * math.sin(math.radians(angle)))
    return intersection_point(wall_coordinates, [cg,arb_point])

#plot-related
def plot_wedge(active_failure_wedge, w_type):
    """
    plots a wedge object
    :param active_failure_wedge: wedge object
    :param w_type: wedge type: 0 if not the failure wedge
    :return: None
    """

    soil_x_coor = []
    soil_y_coor = []

    all_coor = list(active_failure_wedge.coordinates)
    soil = all_coor

    for coor in soil:
        soil_x_coor.append(coor[0])
        soil_y_coor.append(coor[1])

    if w_type == 0: plt.plot(soil_x_coor, soil_y_coor, "grey", linewidth=0.5)
    else: plt.plot(soil_x_coor, soil_y_coor, "brown", linewidth=2.0)


def plot_location(cg, p_active_point, p_active):
    plt.plot([cg[0], p_active_point[0]], [cg[1], p_active_point[1]], "green", linewidth=1, linestyle="dashed")
    plt.plot([cg[0]], [cg[1]], marker='o', markersize=5, color="brown")
    plt.plot([p_active_point[0]], [p_active_point[1]], marker='o', markersize=7, color="deepskyblue")
    plt.plot([p_active_point[0]], [p_active_point[1]], marker='x', markersize=7, color="brown")

    plt.text(p_active_point[0] - 0.5, p_active_point[1], "{:.2f} @ ({:.2f}, {:.2f})".format(p_active, p_active_point[0], p_active_point[1]), ha='right', va='center')

def plot_lines(line_loads, wall_coordinates, friction_angle, failure_angle, backfill_angle, dista, init_soil_coor, h_cr, result):

    it = 0
    p_lines = [(x) for x in result[2:]]
    p_lines.insert(0,0)
    for line in line_loads:
        it+=1
        line_coor = (init_soil_coor[0] + line[1] * math.cos(math.radians(backfill_angle)), h_cr + init_soil_coor[1] + line[1] * math.sin(math.radians(backfill_angle)))

        beta_intersection= [0, 0]
        if line[1] * math.cos(math.radians(backfill_angle)) <= dista:
            beta_intersection = get_line_intersection_with_backfill(wall_coordinates,failure_angle,line_coor)

        plt.plot([line_coor[0], beta_intersection[0]], [line_coor[1], beta_intersection[1]], "purple", linewidth=0.75, linestyle="dashed")

        phi_intersection = get_line_intersection_with_backfill(wall_coordinates,friction_angle,line_coor)
        plt.plot([line_coor[0], phi_intersection[0]], [line_coor[1], phi_intersection[1]], "purple", linewidth=0.75, linestyle="dashed")


        point_of_delta_p = ((beta_intersection[0] + phi_intersection[0]*2)/3, (beta_intersection[1] + phi_intersection[1]*2)/3)

        plt.plot([point_of_delta_p[0]], [point_of_delta_p[1]], marker='o', markersize=5, color="deepskyblue")
        plt.plot([point_of_delta_p[0]], [point_of_delta_p[1]], marker='x', markersize=5, color="brown")

        plt.text(point_of_delta_p[0]-0.5, point_of_delta_p[1], "{:.2f}".format(p_lines[it]-p_lines[it-1])+" @ ({:.2f}, {:.2f})".format(point_of_delta_p[0], point_of_delta_p[1]), ha='right', va='center')






#Analytical Solver
def calculate_p_active(weight, cohesion, adhesion, wedge_angle, wall_angle, int_friction, wall_friction):
    """
    Analytical solution to Columb's Force Polygon
    :param weight: weight of wedge
    :param cohesion: cohesion force on soil slip plane
    :param adhesion: adhesion force on interface
    :return: Pactive
    """
    k = weight - cohesion*math.sin(math.radians(wedge_angle)) - adhesion*math.sin(math.radians(wall_angle))
    m = - cohesion*math.cos(math.radians(wedge_angle)) + adhesion*math.cos(math.radians(wall_angle))
    omega = 90 - wedge_angle + int_friction
    phi = 90 - wall_angle + wall_friction

    return (k + m * math.tan(math.radians(omega))) \
           / (math.cos(math.radians(phi))*math.tan(math.radians(omega))+math.sin(math.radians(phi)))

def get_inputs():
    """gets input from user"""
    global INPUT_MSGS

    inputs = []
    for msg in INPUT_MSGS:
        good_input = False
        while not good_input:
            try:
                inputted = input(msg[1])
                if inputted == "q": quit()
                elif msg[0] == "int": inputs.append(int(inputted))
                elif msg[0] == "float": inputs.append(float(inputted))
                elif msg[0] == "line": inputs.append(format_lines(inputted))
                good_input = True
            except: pass

    return inputs

#Display
def line_loads_display(line_loads):
    final = """"""
    for line in line_loads:
        final += (" "*56 + str(line[0]) + " @ " +str(line[1])+"\n")
    return final
def show_data():
    return """
                You have input:
                00: number of wedges,                   {:.0f}
                01: inclined step width (m),            {:.2f}
                02: wall angle (deg),                   {:.2f}
                03: backfill angle (deg),               {:.2f}
                04: vertical wall height (m),           {:.2f}
                05: density (kN/m3),                    {:.2f}
                06: cohesion (kPa),                     {:.2f}
                07: adhesion (kPa),                     {:.2f}
                08: internal friction angle (deg),      {:.2f}
                09: wall friction angle (deg),          {:.2f}
                10: uniform surcharge (kPa),            {:.2f}
                11: GWT level (m),                      {:.2f}\n""".format(DATA[0], DATA[1], DATA[2], DATA[3], DATA[4],
                                                                       DATA[5], DATA[6], DATA[7], DATA[8], DATA[9], DATA[10], DATA[11]) \
           + " " * 16+"""12: Line Loads (kN/m') @ distance (m), \n"""+line_loads_display(DATA[12])

def solve_and_present():
    global DATA
    solution = main_function(DATA)
    clear()
    warning = """"""
    if solution[0][1] > 90:
        warning += "\n[WARNING] FAILURE ANGLE TOO LARGE\n"
    if solution[0][0] == DATA[0]*DATA[1]:
        warning = warning + "\n" """\n[WARNING] ESTIMATED ACTIVE FAILURE WEDGE IS THE LAST WEDGE! \n"""+ 10*" "+ """YOU SHOULD USE MORE WEDGES!\n"""


    if DATA[12][0][0] == 0 and len(DATA[12])  == 1:
        return warning + """
        *********************************************************
        ************************ SOLUTION ***********************
        *********************************************************
        
        FAILURE WEDGE INCLINED DISTANCE: {:.3f} m
        FAILUER WEDGE ANGLE:             {:.3f} deg
        ACTIVE EARTH PRESSURE:           {:.3f} kN/m'
        
        *********************************************************
        """.format(solution[0][0], solution[0][1], solution[1]) + "\n" * 4 + show_data()
    else:
        lls = ""
        for i in range(len(solution[2:])):
            if i == 0: sol = solution[i+2]
            else: sol = solution[i+2] - solution[i+1]
            lls += " "*49+ "{:.3f}".format(sol)+"\n"


        return warning + """
        *********************************************************
        ************************ SOLUTION ***********************
        *********************************************************
        
        NATURAL FAILURE WEDGE INCLINED DISTANCE: {:.3f} m
        FAILUER WEDGE ANGLE:                     {:.3f} deg
        ACTIVE EARTH PRESSURE:                   {:.3f} kN/m'
        
        DELTA ACTIVE PRESSURES:                  kN/m'\n""".format(solution[0][0], solution[0][1], solution[1]) + lls + """
        
        *********************************************************""" + "\n" * 4 + show_data()




###################################################
#MAIN FUNCTION
def main_function(data, water_density=10):
    """
    UNITS IN m and KN
    :param no_wedges: number of wedges
    :param inclined_width: incremental inclined width on the backfill
    :param wall_angle: horizontal angle (clockwise)
    :param backfill_angle: horizontal angle (anti-clockwise)
    :param vertical_wall_height
    :param density
    :param cohesion: soil cohesion (shear param)
    :param adhesion: soil adhesion to wall (interface param)
    :param int_friction: internal friction angle (shear param)
    :param wall_friction: delta; wall friction angle (interface param)
    :param uniform_surcharge=0,
    :param GWT_level: elev. from base. -1 for no water.
    :param line_load: tuple (magnitude, inclined distance)
    :return: Lateral Earth Pressure in KN
    """

    #Data Definitions
    [no_wedges,
    inclined_width,
    wall_angle,
    backfill_angle,
    vertical_wall_height,
    density,
    cohesion,
    adhesion,
    int_friction,
    wall_friction,
    uniform_surcharge,
    GWT_level,
    line_load] = data


    # define wall points coordinates
    wall_coordinates = [(-1*vertical_wall_height/(math.tan(math.radians(wall_angle))), vertical_wall_height), (0,0)]
    adhesion_force = adhesion * distance(wall_coordinates[0], wall_coordinates[1])


    #hcrack
    ka = (1 - math.sin(math.radians(int_friction))) / (1 + math.sin(math.radians(int_friction))) #Rankine's
    if GWT_level == -1:
        h_crack = (2 * cohesion - uniform_surcharge * math.sqrt(ka)) / (density * math.sqrt(ka))
    else:
        stress_at_GWT = ka * (density * (vertical_wall_height - GWT_level) + uniform_surcharge) - 2 * cohesion * math.sqrt(ka)  # also Rankine's
        if stress_at_GWT >= 0:
            h_crack = (2 * cohesion - uniform_surcharge * math.sqrt(ka)) / (density * math.sqrt(ka))
        else:
            h_crack = -1 * stress_at_GWT / (ka * (density - water_density))

    #abort if h_crack > vertical height
    if h_crack > vertical_wall_height:
        print("h_crack is larger than the vertical wall height! recheck your input..")
        return -1

    if h_crack < 0: h_crack = 0

    #init coordinates
    init_inclination = h_crack*(math.sin(math.radians(90-wall_angle))/math.sin(math.radians(wall_angle + backfill_angle)))
    init_soil_coor = (wall_coordinates[0][0] + init_inclination * math.cos(math.radians(backfill_angle)),
                      wall_coordinates[0][1] + init_inclination * math.sin(math.radians(backfill_angle))-h_crack)

    p_active = float('-inf')
    active_failure_wedge = None
    active_failure_angle = 0
    dista = 0
    for num in range(1, no_wedges+1):
        iter_coor = (init_soil_coor[0]+ num*inclined_width*math.cos(math.radians(backfill_angle)),
                     init_soil_coor[1]+ num*inclined_width*math.sin(math.radians(backfill_angle)))

        failure_angle = line_angle(wall_coordinates[1], iter_coor)
        if failure_angle < 0: failure_angle = failure_angle + 180

        top_iter_coor = (iter_coor[0], iter_coor[1] + h_crack)

        #determine coordinates
        if h_crack > 0: coordinates = [wall_coordinates[0], init_soil_coor, wall_coordinates[1], iter_coor, top_iter_coor]
        else: coordinates = [wall_coordinates[0], wall_coordinates[1], iter_coor]


        wedge = FailureWedge(density, coordinates, GWT_level, water_density)
        cohesion_force = cohesion * distance(iter_coor, wall_coordinates[1])
        surcharge_force = uniform_surcharge * num * inclined_width
        p_active_current = calculate_p_active(wedge.weight+surcharge_force, cohesion_force, adhesion_force, failure_angle, wall_angle, int_friction, wall_friction)

        plot_wedge(wedge, 0)

        if p_active_current > p_active:
            p_active = p_active_current
            active_failure_wedge = wedge
            active_failure_angle = failure_angle
            dista = num * inclined_width


    cg_wedge = active_failure_wedge.get_wedge_cg(uniform_surcharge, dista, init_soil_coor, backfill_angle)
    p_active_point = get_p_active_point(wall_coordinates,active_failure_angle,cg_wedge)

    plt.plot([init_soil_coor[0], iter_coor[0]], [init_soil_coor[1], iter_coor[1]], "grey", linewidth=0.75, linestyle="dashed")
    plt.plot([wall_coordinates[0][0], iter_coor[0]], [wall_coordinates[0][1], iter_coor[1]+h_crack], "saddlebrown", linewidth=2)

    plot_wedge(active_failure_wedge, 1)




    GWT_intersection = intersection_point(wall_coordinates, [(0, GWT_level), (100, GWT_level)])
    if GWT_level != -1: plt.plot([GWT_intersection[0], dista], [GWT_level, GWT_level], "blue", linewidth=1)

    if surcharge_force != 0:
        plt.plot([wall_coordinates[0][0], wall_coordinates[0][0], top_iter_coor[0], top_iter_coor[0]], [wall_coordinates[0][1], wall_coordinates[0][1]+0.2, top_iter_coor[1]+0.2, top_iter_coor[1]], "purple", linestyle="dotted", linewidth=1.5)
        plt.text(init_soil_coor[0]+ dista*math.cos(math.radians(backfill_angle)), init_soil_coor[1]+ dista*math.sin(math.radians(backfill_angle)) + h_crack + 1, str(uniform_surcharge) + " kN/m2", ha='left')


    plt.plot([wall_coordinates[0][0], wall_coordinates[1][0], wall_coordinates[0][0]],
             [wall_coordinates[0][1], wall_coordinates[1][1], 0], "black", linewidth=2.5)

    plot_location(cg_wedge, p_active_point, p_active)
    plt.text(0, -0.25, "(0, 0)", ha='center', va='top')

    result = [[dista, active_failure_angle], p_active]

    if (line_load[0][0] == 0) and (len(line_load) == 1) : return result

    else:
        for i in range(len(line_load)):
            p_plus = float('-inf')
            for num in range(1, no_wedges + 1):
                iter_coor = (init_soil_coor[0] + num * inclined_width * math.cos(math.radians(backfill_angle)),
                             init_soil_coor[1] + num * inclined_width * math.sin(math.radians(backfill_angle)))

                failure_angle = line_angle(wall_coordinates[1], iter_coor)
                if failure_angle < 0: failure_angle = failure_angle + 180

                # determine coordinates
                if h_crack > 0:
                    top_iter_coor = (iter_coor[0], iter_coor[1] + h_crack)
                    coordinates = [wall_coordinates[0], init_soil_coor, wall_coordinates[1], iter_coor, top_iter_coor]

                else:
                    coordinates = [wall_coordinates[0], wall_coordinates[1], iter_coor]

                wedge = FailureWedge(density, coordinates, GWT_level, water_density)
                cohesion_force = cohesion * distance(iter_coor, wall_coordinates[1])
                surcharge_force = uniform_surcharge * num * inclined_width
                total = wedge.weight + surcharge_force

                for j in range(i+1):
                    if line_load[j][1] <= num*inclined_width: total += line_load[j][0]

                p_active_current = calculate_p_active(total, cohesion_force, adhesion_force,
                                                      failure_angle, wall_angle, int_friction, wall_friction)

                if p_active_current > p_plus:
                    p_plus = p_active_current

            delta = p_plus - result[1]

            point_load_coor = (init_soil_coor[0] + line_load[i][1] * math.cos(math.radians(backfill_angle)),
                             init_soil_coor[1] + line_load[i][1] * math.sin(math.radians(backfill_angle))+h_crack)

            plt.text(point_load_coor[0], point_load_coor[1]+1, str(line_load[i][0])+" kN/m'",  ha='center')

            plt.plot([point_load_coor[0], point_load_coor[0]],[point_load_coor[1]+1,point_load_coor[1]], "green", linestyle="dotted")

            result.append(delta)

        plot_lines(line_load, wall_coordinates, int_friction, active_failure_angle, backfill_angle, dista, init_soil_coor, h_crack, result)

        return result


###################################################

#initializer function
def initializer():
    global DATA

    clear()
    print("""
***********************************************************************************************************
************************** WELCOME TO COLUMB ANALYTICAL 1.3.3 | Abdalla Talaat\xa9 ***************************
***********************************************************************************************************
           \n\n
           """)

    print(show_data())

    a = input("input (n) to add new data. input (q) to quit. Press Enter to solve using existing data.  ")
    done = False
    while not done:
        if a == "n":
            try:
                DATA = list(get_inputs())
                done = True
            except:
                quit()

        elif a == "q": quit()
        else:
            done = True
            print(solve_and_present())

#Main loop
while(True):

    initializer()
    answered = False
    while not answered:
        q = input("\nWould you like to calculate again (y/n)?\n Or would you like to edit (e)?\n Or would you like to plot the result (p)?  ")

        if q == "y":
            answered = True
            pass
        elif q == "n":
            clear()
            quit()
        elif q=="e":
            clear()
            print(show_data())
            try:
                edi = int(input("""
                Choose what to edit:
                0: number of wedges,
                1: inclined step width,
                2: wall angle,
                3: backfill angle,
                4: vertical wall height,
                5: density,
                6: cohesion,
                7: adhesion,
                8: internal friction angle,
                9: wall friction angle,
                10: uniform surcharge,
                11: GWT level
                12: Line Loads
                """))

                if edi == 0:
                    DATA[0] = int(input("Input the new value:  "))
                elif edi == 12:
                    input_line_loads = [float(x) for x in
                                        input(
                                            "\nselect consecutively the desired line load(s) and inclined distance(s)\n"
                                            "from R.Wall top separated by a spaces (0 0 for no line load) (kN/m'): ").split()]
                    line_loads = []
                    for i in range(0, len(input_line_loads), 2):
                        line_loads.append(tuple([input_line_loads[i], input_line_loads[i + 1]]))
                    DATA[12] = line_loads
                else:
                    DATA[edi] = float(input("Input the new value:  "))
            except:
                print("BAD INPUT")


        elif q=="p":
            plt.close()

            plt.figure(figsize=(15, 7), tight_layout=True)

            plt.axis('off')
            plt.gca().set_aspect('equal', adjustable='box')
            txt = solve_and_present()
            print(txt)

            plt.figure(figsize=(4, 5.2), tight_layout=True)
            plt.text(0, 0, txt, ha='right', va='baseline', fontsize='x-small')
            plt.axis('off')


            plt.show()
        else:
            clear()
            print(show_data())
