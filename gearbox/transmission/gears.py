from math import pi, atan, tan, radians, cos, sin, degrees, sqrt, floor, ceil

from gearbox.libs.maths import involute, arcinvolute


class Tool(object):
    """

    :param ha_p: addendum of basic rack of cylindrical gears. Normally: ha_p=m.  [mm]
    :param hf_p: dedendum of basic rack of cylindrical gears. Normally: hf_p=1.2*m. [mm}
    :param rho_fp: root fillet radius of the basic rack for cylindrical gears. [mm]
    :param x: profile shift coefficient
    :param c: the tooth tip/root clearance [mm]
    :param alpha_pn: normal pressure angle of the basic rack for cylindrical gears [°]
    :param nc: Only if calculations according to AGMA are needed
    :param rho_ao: Only if calculations according to AGMA are needed
    :param delta_ao: Only if calculations according to AGMA are needed
    """

    def __init__(self, ha_p, hf_p, rho_fp, c, x, alpha_pn, nc=None, rho_ao=None, delta_ao=None):
        self.ha_p = ha_p
        self.hf_p = hf_p
        self.rho_fp = rho_fp
        self.c = c
        self.x = x
        self.alpha_pn = alpha_pn
        self.nc = nc
        self.rho_ao = rho_ao
        self.delta_ao = delta_ao
        


class Material(object):
    """

    :param sigma_h_lim: Allowable hertzian stress [MPa]
    :param sigma_f_lim: Allowable bending stress [MPa]
    :param brinell: Brinell hardness of the tooth flanks.
    :param classification:
            'St'            Wrought normalized low carbon steels
            'St(cast)'      Cast steels
            'GTS(perl)'     Black malleable cast iron (perlitic structure)
            'GGG(ferr)'     Nodular cast iron (perlitic, Cast iron materials bainitic, ferritic structure)
            'GG'            Grey cast iron
            'V'             Through-hardened wrought steels Carbon steels, alloy steels
            'V(cast)'       Through-hardened cast steels Carbon steels, alloy steels
            'Eh'            Case-hardened wrought steels
            'IF'            Flame or induction hardened wrought or cast steels
            'NT(nitr)'      Nitrided wrought steels - Nitriding steels
            'NV(nitr)'      Nitrided wrought steels - Through-hardening steels
            'NV(nitrocar)'  Wrought steels, nitrocarburized Through hardening steels
            
    :param name: Material designation
    :param e: Elastic modulus [MPa]
    :param poisson: Poisson ratio of material
    :param density: Material density [kg/mm**3]
    """

    def __init__(self, sigma_h_lim, sigma_f_lim, brinell, classification, name='', e=206000., poisson=0.3, density=7.83e-6):
        self.sigma_h_lim = sigma_h_lim
        self.sigma_f_lim = sigma_f_lim
        self.classification = classification
        self.name = name
        self.e = e
        self.poisson = poisson
        self.density = density
        self.brinell = brinell


class Lubricant(object):
    """

    :param v40: Viscosity at 40°C, Unit: mm**2/s = cSt
    :param name: Designation
    """

    def __init__(self, v40, name=''):
        self.name = name
        self.v40 = v40


class Gear(object):
    """

    :param profile:
    :param material:
    :param z:
    :param beta:
    :param b:
    :param bs:
    :param alpha:
    :param m:
    :param x:
    :param sr:
    :param rz:
    :param precision_grade:
    :param shaft_diameter:
    :param schema:
    :param l:
    :param s:
    :param backlash:
    :param gear_crown:
    :param gear_condition:
    :param helix_modification:
    """

    def __init__(self, **kw):
        if 'gear' in kw:
            self.profile = kw['gear'].profile
            self.material = kw['gear'].material
            self.z = kw['gear'].z
            self.beta = kw['gear'].beta
            self.alpha = kw['gear'].alpha
            self.m = kw['gear'].m
            self.x = kw['gear'].x
            self.b = kw['gear'].b
            self.bs = kw['gear'].bs
            self.sr = kw['gear'].sr
            self.rz = kw['gear'].rz
            self.precision_grade = kw['gear'].precision_grade
            self.shaft_diameter = kw['gear'].shaft_diameter
            self.schema = kw['gear'].schema
            self.l = kw['gear'].l
            self.s = kw['gear'].s
            self.backlash = kw['gear'].backlash

            self.gear_crown = kw['gear'].gear_crown
            self.helix_modification = kw['gear'].helix_modification

            self.gear_condition = kw['gear'].gear_condition

        else:
            self.profile = kw['profile']
            self.material = kw['material']
            self.z = kw['z']
            self.beta = kw['beta']
            self.alpha = kw['alpha']
            self.m = kw['m']
            self.x = kw['x']
            self.b = kw['b']
            self.bs = kw['bs']
            self.sr = kw['sr']
            self.rz = kw['rz']
            self.precision_grade = kw['precision_grade']
            self.shaft_diameter = kw['shaft_diameter']
            self.schema = kw['schema']
            self.l = kw['l']
            self.s = kw['s']
            self.backlash = kw['backlash']

            self.gear_crown = kw['gear_crown']
            self.helix_modification = kw['helix_modification']

            self.gear_condition = kw['gear_condition']

        self.xmin = (1 - sqrt(self.z * sin(radians(self.alpha)))) / 2
        self.alpha_t = degrees(atan(tan(radians(self.alpha)) / cos(radians(self.beta))))
        self.addendum = self.m * (1 + self.x)
        self.dedendum = self.m * (1 - self.x) + self.profile.c
        self.h = self.dedendum + self.addendum
        self.d = self.m * self.z / cos(radians(self.beta))
        self.da = self.d + 2*self.addendum   
        self.df = self.d - 2*self.dedendum
        self.db = self.d * cos(radians(self.alpha_t))
        self.rho_f = self.profile.rho_fp * self.m
        self.mt = self.m / cos(radians(self.beta))
        self.p_b = self.m * cos(radians(self.alpha)) * pi
        self.p_n = pi * cos(radians(self.alpha))

        # FIXME self.sn = ((pi/2) + 2 * 0 *self.m*self.x*tan(self.Alpha))
        self.beta_b = degrees(atan(self.db * tan(radians(self.beta)) / self.d))

        self.zn = self.z / (cos(radians(self.beta)) * cos(radians(self.beta_b)) ** 2)

        mc = self.__interval_calc([0, 0.5, 2.0, 3.5, 6.0, 25, 40, 70], self.m)
        dc = self.__interval_calc([0, 5, 20, 50, 125, 280, 560, 1000, 1600, 2500, 4000, 6000, 8000, 10000], self.d)
        bc = self.__interval_calc([0, 4, 10, 20, 40, 80, 160, 250, 400, 650, 1000], self.b)

        f_pt = 0.3 * (mc + 0.4 * sqrt(dc)) + 4.
        f_p = 0.3 * mc + 1.25 * sqrt(dc) + 7.
        f_a = 3.2 * sqrt(mc) + 0.22 * sqrt(dc) + 0.27
        f_beta = 0.1 * sqrt(dc) + 0.63 * sqrt(bc) + 4.2
        f_f_alpha = 2.5 * sqrt(mc) + 0.17 * sqrt(dc) + 0.5
        f_h_alpha = 2 * sqrt(mc) + 0.14 * sqrt(dc) + 0.5
        f_h_beta = 0.07 * sqrt(dc) + 0.45 * sqrt(bc) + 3.

        self.f_pt = self.__q(f_pt, self.precision_grade)
        self.f_p = self.__q(f_p, self.precision_grade)
        self.f_alpha = self.__q(f_a, self.precision_grade)
        self.f_beta = self.__q(f_beta, self.precision_grade)
        self.f_f_alpha = self.__q(f_f_alpha, self.precision_grade)
        self.f_h_alpha = self.__q(f_h_alpha, self.precision_grade)
        self.f_h_beta = self.__q(f_h_beta, self.precision_grade)
        self.f_f_beta = self.__q(f_h_beta, self.precision_grade)
        self.f_h_beta6 = self.__q(f_h_beta, 6)

    @staticmethod
    def __interval_calc(intervals, attr):
        for i in range(1, intervals.__len__()):
            if intervals[i - 1] <= attr < intervals[i]:
                return sqrt(intervals[i] * intervals[i - 1])

    @staticmethod
    def __q(x, qiso):
        x *= 2 ** (0.5 * (qiso - 5))
        if x >= 10:
            x = round(x)
        elif 5 <= x < 10:
            if x % 1 <= 0.25 or (0.5 <= x % 1 % 1 <= 0.75):
                x = floor(x * 2) * 0.5
            else:
                x = ceil(x * 2) * 0.5
        else:
            x = round(x, 1)
        return x


class Transmission(object):
    """
    
    :param lubricant:
    :param rpm_in: the rotation speed of the first gear
    :param gear_box_type: used for AGMA calculation
    :param p: the transmitted power in kW
    :param l: the design life in hours
    :param gears: tuple of Gear objects
    :param ka: the application factor
    :param sf_min: the minimum safety factor for tooth bending
    :param sh_min: the minumum safety factor for flank pitting
    :param favorable_contact: Can be set True if the contact pattern is known to be good.
    """

    def __init__(self, **kw):

        gears = kw['gears']
        self.rpm_in = kw['rpm_in']
        self.rpm_out = gears[0].z/gears[1].z*self.rpm_in
        self.ka = kw['ka']
        self.sh_min = kw['sh_min']
        self.sf_min = kw['sf_min']

        self.v40 = kw['lubricant'].v40
        self.gear_box_type = kw['gear_box_type']

        self.u = self.rpm_out / self.rpm_in
        self.p = kw['p']
        self.l = kw['l']
        self.favorable_contact = kw['favorable_contact']

        if 'k_h_alpha' in kw:
            self.k_h_alpha = kw['k_h_alpha']
        else:
            self.k_h_alpha = None

        if 'k_h_beta' in kw:
            self.k_h_beta = kw['k_h_beta']
        else:
            self.k_h_beta = None
        
        if 'k_f_beta' in kw:
            self.k_f_beta = kw['k_f_beta']    
        else:
            self.k_f_beta = None
            
        self.pair = self.__calculate(gears[0], gears[1])

    def __calculate(self, gear_one, gear_two):
        if gear_one.m is not gear_two.m:
            raise Exception("the modulus of the two gears most be equal")
        else:
            self.m = gear_one.m

        if gear_one.alpha is not gear_two.alpha:
            raise Exception("the pressure angle of the two gears most be equal")
        else:
            self.alpha = gear_one.alpha
            self.alpha_t = gear_one.alpha_t

        self.u = gear_two.z / gear_one.z
        inv = involute(gear_one.alpha_t) + 2 * (gear_one.x + gear_two.x) / (gear_one.z + gear_two.z) * tan(
            radians(gear_one.alpha))
        self.alpha_wt = arcinvolute(inv)
        self.a = ((gear_one.z + gear_two.z) * gear_one.m) / (2 * cos(radians(gear_one.beta)))
        self.aw = self.a * cos(radians(gear_one.alpha)) / cos(radians(self.alpha_wt))
        self.epsilon_alpha = (
            0.5 * (sqrt(gear_one.da ** 2 - gear_one.db ** 2) + sqrt(gear_two.da ** 2 - gear_two.db ** 2)) - self.a * sin(
                radians(self.alpha_wt))) / (pi * gear_one.m * cos(radians(gear_one.alpha_t)) / cos(radians(gear_one.beta)))
        self.epsilon_beta = gear_one.b * sin(radians(gear_one.beta)) / (gear_one.m * pi)
        self.epsilon_gama = self.epsilon_alpha + self.epsilon_beta
        self.v = self.rpm_in * gear_one.d * pi / 60000
        self.ft = 1000. * self.p / self.v
        if self.ka * self.ft / gear_one.b < 100:
            self.fmt = 100
        else:
            self.fmt = self.ka * self.ft / gear_one.b

        self.xsum = gear_one.x + gear_two.x
        self.gear_one = gear_one
        self.gear_two = gear_two
