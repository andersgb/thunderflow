import thulip_ext as th
import numpy as np
import logging
logging.basicConfig(level=logging.ERROR, format='%(asctime)s - %(levelname)s - %(message)s')


class pythulip:
    def __init__(self, model_name, components):
        self.t = th.create_thulip(model_name, components)

    def __getitem__(self, key):
        return th.thermo_get(self.t, key)
    
    def __setitem__(self, key, value):
        if isinstance(value, list):
            th.set_by_pylist(self.t, key, value);
        else:
            self.t.set_scalar(key, value, 0)
    
    def calc(self):
        self.t.calc();

def set_SI_units(mod):
    mod["unit_1A"] = mod["si_e_current"]
    mod["unit_1dm"] = mod["si_length"]
    mod["unit_1e3k"] = mod["si_temperature"]
    mod["unit_1e4ton"] = mod["si_mass"]
    mod["unit_1mol"] = mod["si_amount"]
    mod["unit_1s"] = mod["si_time"]

def iterate_xy(mod, var_a, value_a, var_b, value_b):
    """ Newton iterate on X and Y to obtain value_a and value_b"""
    while True:
        jac = np.array([[ mod[var_a+'_x'][0], mod[var_a+'_y'][0] ],
                        [ mod[var_b+'_x'][0], mod[var_b+'_y'][0] ]])
        b = np.array([ value_a - mod[var_a][0], value_b - mod[var_b][0]])
        x = np.linalg.solve(jac, b)
        
        logging.debug(str(x))
        mod['delta_x'] = x[0]
        mod['delta_y'] = x[1]
        mod['delta_z'] = [0, 0]
        #TODO: step length control
        mod['delta_step'] = 1.0
        mod.calc()
        mod['var_x'] = mod['state_x']
        mod['var_y'] = mod['state_y']

        logging.debug( str(mod['state_t']) )
        logging.debug( str(mod['state_p']) )
        logging.debug( str(mod['delta_norm']) )
        logging.debug('--')
        if mod['delta_norm'][0] < 1e-10:
            mod['delta_step'] = 0.0
            break

def iterate_x(mod, var_a, value_a):
    mod['delta_y'] = 0
    mod['delta_z'] = [0.0 for _ in xrange(mod.t.number_of_elements('delta_z'))]
    while True:
        jac = mod[var_a+'_x'][0]
        b = value_a - mod[var_a][0]
        mod['delta_x'] = b/jac
        mod['delta_step'] = 1.0
        mod.calc()
        mod['var_x'] = mod['state_x']
        
        if mod['delta_norm'][0] < 1e-10:
            mod['delta_step'] = 0.0
            break
        
### TEEESSST
# p = pythulip('pr_c4', ['methane', 'ethane'])
# p.calc()
# p['var_n'] = [0.5, 0.5]
# p['var_v'] = p['fix_rgas'][0]*0.300/0.001;
# p.calc()
# iterate_xy(p, 'state_t', 0.300, 'state_p', 0.001)
