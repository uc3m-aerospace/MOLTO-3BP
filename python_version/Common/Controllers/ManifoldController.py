import numpy as np
# TODO Skip following import
from Utility.PCR3BP_helpers import PCR3BP_propagator, PCR3BP_state_derivs
from Utility.manifold_helpers import construct, plotm, fourierTest, load_variables, query_func, query_func_text
import Utility.helpers as h
from Common.Controllers.Halo.HaloController import HaloController


class ManifoldController:
    @staticmethod
    def main_runner(data):
        h.log('\nManifolds Propagator Software\n')

        ## Initialize variables
        data = load_variables(data)

        ## 1. Orbit family
        if data['input_seq']:
            if not data['text']:
                Sequence = query_func()
                h.log('\n')
                h.log(Sequence)
            else:
                Sequence = query_func_text(data)
                h.log('\n')
                h.log(Sequence)


            data['flags'] = [1, 1, 1]
            data['tf'] = 4
            data['nmax'] = 50
            data['tol'] = 1e-15

            orbitDef = []

            for i in range(1, Sequence['it'] + 1):
                h.log('\nEvaluating Orbit ' + str(i))
                data['type'] = Sequence['type' + str(i)]

                if data['type'] == 'LY':
                    data['Ax_tgt'] = Sequence['Ax' + str(i)]
                    data['LP'] = Sequence['LP' + str(i)]
                    # Todo Compute Here for Lyupanov
                    # ly_controller = LyupanovController(data)
                    # orbitDef.append(ly_controller.main())
                    # orbitDef.append(Lyapunov(data))
                else:
                    data['Az'] = Sequence['Az' + str(i)]
                    data['phi'] = Sequence['phi' + str(i)]
                    data['m'] = Sequence['m' + str(i)]
                    data['LP'] = Sequence['LP' + str(i)]
                    halo_controller = HaloController(data)
                    orbitDef.append(halo_controller.main())

        else:
            if data['type'] == 'LY':
                from Common.Controllers.Lyapunov.LyapunovController import Lyapunov
                orbitDef = Lyapunov(data)

            elif data['type'] == 'HL':
                data['flags'] = [1, 1, 1]
                data['tf'] = 4
                data['nmax'] = 50
                data['tol'] = 1e-15
                from Halo_Orbits.Halo_Main import Halo_Main
                orbitDef = Halo_Main(data)

            else:
                raise Exception('Manifolds_Main:typeError.' + \
                                '    The type selected is not valid [\'LY\'][\'HL\']!')

        ## 2.1 Construct manifolds
        npoints = data['npoints']  # Number of iterations = npoints*2

        stop_fun = h.poinc_crossing

        if data['input_seq']:
            for j in range(1, Sequence['it'] + 1):

                states_s = []
                SF_s = np.array([]).reshape(-1, 4)
                states_u = []
                SF_u = np.array([]).reshape(-1, 4)

                if 'ang' + str(j) not in Sequence:
                    break

                ang = Sequence['ang' + str(j)] % 360

                if Sequence['LP' + str(j)] - 1:
                    angmin = min(np.arctan2(orbitDef[j - 1][0][1], orbitDef[j - 1][0][0] \
                                            - data['params'][0]))
                    if ang * np.pi / 180 > angmin + 2 * np.pi or ang * np.pi / 180 < -angmin:
                        h.log('This angle is not valid (intersecting initial orbit)!')
                        if abs(ang * np.pi / 180 - (angmin + 2 * np.pi)) < abs(ang * np.pi / 180 + angmin):
                            ang = 4 / 3 * angmin * 180 / np.pi + 360
                            h.log('Reassigning to %3.2f...' % ang)
                        else:
                            ang = -4 / 3 * angmin * 180 / np.pi
                            h.log('Reassigning to %3.2f...' % ang)
                else:
                    angmin = min(np.arctan2(orbitDef[j - 1][0][1], data['params'][0] \
                                            - orbitDef[j - 1][0][0])) \
                             + np.pi
                    if ang * np.pi / 180 > angmin and ang * np.pi / 180 < -angmin + 2 * np.pi:
                        h.log('This angle is not valid (intersecting initial orbit)!')
                        if abs(ang * np.pi / 180 - (-angmin + 2 * np.pi)) < abs(ang * np.pi / 180 - angmin):
                            ang = -4 / 3 * angmin * 180 / np.pi + 420
                            h.log('Reassigning to %3.2f...' % ang)
                        else:
                            ang = 4 / 3 * angmin * 180 / np.pi - 60
                            h.log('Reassigning to %3.2f...' % ang)

                h.log('\nConstructing Manifolds...\n')
                h.log('Poincaré section angle = %3.1f' % ang)

                h.log('Unstable manifolds...')
                [states_s_e, times_s_e, SF_s_e, states_u, times_u, SF_u] = construct(
                    data['params'], orbitDef[j - 1], data['prnt_out_dt'], npoints, 1,
                    -1, stop_fun, ang, data['pos'][Sequence['LP' + str(j)] - 1][0])

                if 'Ax' + str(j + 1) in Sequence or 'Az' + str(j + 1) in Sequence:
                    h.log('Stable manifolds...')
                    [states_s, times_s, SF_s, states_u_e, times_u_e, SF_u_e] = construct(
                        data['params'], orbitDef[j], data['prnt_out_dt'], npoints, 1,
                        1, stop_fun, ang, data['pos'][Sequence['LP' + str(j + 1)] - 1][0])

                    h.log('\nPlotting manifolds\n')

                    ## 2.3 Plot manifolds

                    pos = np.append(data['pos'][Sequence['LP' + str(j)] - 1],
                                    data['pos'][Sequence['LP' + str(j + 1)] - 1])

                    plotm(data['params'][0], data['params'][1], pos, orbitDef[j - 1:j + 1],
                          states_s, SF_s, states_u, SF_u, ang, angmin)

                else:

                    h.log('\nPlotting manifolds\n')

                    ## 2.3 Plot manifolds
                    plotm(data['params'][0], data['params'][1],
                          data['pos'][Sequence['LP' + str(j)] - 1], orbitDef[j - 1],
                          states_s, SF_s, states_u, SF_u, ang, angmin)

        else:
            for i in data['poincSec']:

                ang = i % 360
                if data['LP'] - 1:
                    angmin = min(np.arctan2(orbitDef[0][1], orbitDef[0][0] - data['params'][0]))
                    if ang * np.pi / 180 > angmin + 2 * np.pi or ang * np.pi / 180 < -angmin:
                        h.log('This angle is not valid (intersecting initial orbit)!')
                        if abs(ang * np.pi / 180 - (angmin + 2 * np.pi)) < abs(ang * np.pi / 180 + angmin):
                            ang = 4 / 3 * angmin * 180 / np.pi + 360
                            h.log('Reassigning to %3.2f...' % ang)
                        else:
                            ang = -4 / 3 * angmin * 180 / np.pi
                            h.log('Reassigning to %3.2f...' % ang)
                else:
                    angmin = min(np.arctan2(orbitDef[0][1], data['params'][0] - orbitDef[0][0])) \
                             + np.pi
                    if ang * np.pi / 180 > angmin and ang * np.pi / 180 < -angmin + 2 * np.pi:
                        h.log('This angle is not valid (intersecting initial orbit)!')
                        if abs(ang * np.pi / 180 - (-angmin + 2 * np.pi)) < abs(ang * np.pi / 180 - angmin):
                            ang = -4 / 3 * angmin * 180 / np.pi + 420
                            h.log('Reassigning to %3.2f...' % ang)
                        else:
                            ang = 4 / 3 * angmin * 180 / np.pi - 60
                            h.log('Reassigning to %3.2f...' % ang)

                h.log('\nConstructing Manifolds...\n')
                h.log('Poincaré section angle = %3.1f' % ang)

                [states_s, times_s, SF_s, states_u, times_u, SF_u] = construct(
                    data['params'], orbitDef, data['prnt_out_dt'], npoints, data['d'],
                    data['branch'], stop_fun, ang, data['pos'][data['LP'] - 1][0])

                h.log('\nPost-processing data...\n')

                ## 2.2 Fourier analysis
                fourierTest(data['params'][0], data['params'][1], data['pos'][data['LP'] - 1],
                            states_s, states_u, ang, data)

                h.log('\nPlotting manifolds\n')

                ## 2.3 Plot manifolds
                plotm(data['params'][0], data['params'][1], data['pos'][data['LP'] - 1],
                      orbitDef, states_s, SF_s, states_u, SF_u, ang, angmin)
