
'''

'''
import re, os
import numpy as np
import copy

def kenoOut(outName):
    '''
    Extract KENO output, keff, sigma, runtime
    '''
    runtim_pattern        = r'Sequence total time\s+\d+.\d+e*\+\d+'
    keff_pattern          = r'best estimate system k-eff \s+ \d*.\d*\s+\+ or \-\s*\d*.\d*'
    system_total_pattern  = r'system total =\s+\d*.*\d+E*\-*\+*\d*'
    nubar_pattern         = r'system nu bar \s+ \d*.\d*E\+\d*\s+\+ or \-\s*\d*.\d*E\-\d*'
    mfp_pattern           = r'system mean free path \(cm\) \s+ \d*.\d*E\+*\-*\d*\s+\+ or \-\s*\d*.\d*E\-\d*'
    energy                = r'Energy of average lethargy of Fission (eV) k-eff \s+ \d*.\d*\s+\+ or \-\s*\d*.\d*'

    keff                = None, None
    runTime, mfp, nubar = None, None, None
    total_absorption, total_fission, total_leakage = None

    s             = open(outName).read()
    match_keff    = re.findall(keff_pattern, s)
    match_runtime = re.findall(runtim_pattern, s)
    system_total  = re.search(system_total_pattern, s)
    nubar         = re.findall(nubar_pattern, s)
    mfp           = re.findall(mfp_pattern, s)

    if match_keff:
        keff = [
                float(match_keff[-1].split()[4]),
                float(match_keff[-1].split()[-1])
               ]

    if match_runtime:
        timeData = match_runtime[-1].split()[-1]
        runTime  = float(timeData)/60

    if system_total:
        # note: standard deviation is given as percentage so it's converted to absolute
        total_fission    = float(system_total[0].split()[3]), \
                           float(system_total[0].split()[4])/100
        total_absorption = float(system_total[0].split()[5]), \
                           float(system_total[0].split()[6])/100
        total_leakage    = float(system_total[0].split()[7]), \
                           float(system_total[0].split()[8])/100

    if nubar:
        nubar = float(nubar[0].split()[3]), float(nubar[0].split()[7])

    if mfp:
        mfp   = float(mfp[0].split()[5]), float(mfp[0].split()[9])

    return keff[0], keff[1], runTime, mfp, nubar, total_absorption,\
           total_fission, total_leakage


def getKENOReactinData(name):
    '''
    Read KENO .out file to extract system absorptions and fissions 
    '''
    data  = {}
    lines = open(name).readlines()
    for no, line in enumerate(lines):
        if 'system total' in line:
            break

    for line in lines[no+1:]:
        if 'elapsed time' in line:
            break
        elif len(line) > 1 :
            info = line.split()
            if len(info) == 6:
                # note: standard deviation is given as percentage so it's converted to absolute
                unit, total_fissions, total_absorption = info[0],\
                                                         np.array([float(info[2]), float(info[3])/100]),\
                                                         np.array([float(info[4]), float(info[5])/100])
                if unit not in data:
                    data[unit] = {'total_fissions':np.zeros(2),
                                  'total_absorption':np.zeros(2)
                                 }

            elif len(info) == 5:
                total_fissions, total_absorption = (float(info[1]), float(info[2])/100),\
                                                   (float(info[3]), float(info[4])/100)

            data[unit]['total_fissions'] += total_fissions
            # for absorption in unit 1, the pebbel, consider only the fuel region
            if unit == '1':
                if total_fissions[0] == 0:
                    continue
                    #data[unit]['total_absorption'] += np.array([0, 0])
                else:
                    data[unit]['total_absorption'] += total_absorption
            else:
                data[unit]['total_absorption'] += total_absorption

    return data

def getKENOGroupReactinData(name):
    '''
    read KENO .out file to extract system absorptions and fissions 
    '''
    data_in = {
               'fuel' :{'fission': np.zeros(2), 'absorption': np.zeros(2)},
               'total':{'fission': np.zeros(2), 'absorption': np.zeros(2),
                        'leakage': np.zeros(2)}
              }

    data = {
            'thermal': copy.deepcopy(data_in),
            'fast'   : copy.deepcopy(data_in)
           }
    lines = open(name).readlines()

    for no, line in enumerate(lines):
        if 'group  fission  unit   region' in line:

            info = lines[no+3].split()
            group, fission_fraction, total_fissions, total_absorption, leakage = int(info[0]), info[1], \
                                                    np.array([float(info[2]), (float(info[3])/100)**2]),\
                                                    np.array([float(info[4]), (float(info[5])/100)**2]),\
                                                    np.array([float(info[6]), (float(info[7])/100)**2])
                                                    
            if group <= 214:
                g = 'fast'
            elif group > 214:
                g = 'thermal'

            data[g]['total']['absorption'] += total_absorption
            data[g]['total']['fission']    += total_fissions
            data[g]['total']['leakage']    += leakage

            for i in lines[no+4:]:
                if '-  -  -  -  -  -  -  -  -  -  -' in i:
                    break
                info = i.split()
                if len(info) == 6:
                    # note: standard deviation is given as percentage so it's converted to absolute
                    unit, total_fissions, total_absorption = info[0],\
                                                            np.array([float(info[2]), (float(info[3])/100)**2]),\
                                                            np.array([float(info[4]), (float(info[5])/100)**2]),
                    # if unit not in data:
                    #     data[group][unit] = {'total_fissions':np.zeros(2), 'total_absorption':np.zeros(2)}

                elif len(info) == 5:
                    total_fissions, total_absorption = (float(info[1]), (float(info[2])/100)**2),\
                                                       (float(info[3]), (float(info[4])/100)**2)

                #data[group][unit]['total_fissions'] += total_fissions
                # for absorption in unit 1, the pebbel, consider only the fuel region
                if unit == '1':
                    if total_fissions[0] == 0:
                        continue
                        #data[unit]['total_absorption'] += np.array([0, 0])
                    else:
                        data[g]['fuel']['absorption'] += total_absorption
                        data[g]['fuel']['fission'] += total_fissions

            for g in ['thermal', 'fast']:
                for region in ['total', 'fuel']:
                    data[g][region]['absorption'][1] = np.sqrt(data[g][region]['absorption'][1])
                    data[g][region]['fission'][1] = np.sqrt(data[g][region]['fission'][1])

            data['thermal']['total']['leakage'][1] = np.sqrt(data['thermal']['total']['leakage'][1])
            data['fast']['total']['leakage'][1] = np.sqrt(data['fast']['total']['leakage'][1])

    return data


def readKenoFlux(fname):
    
    '''Read Keno flux'''
    flux       = []
    flux_sigma = []
    lines      = open(fname).readlines()

    for no, line in enumerate(lines):
        if 'fluxes for this problem' in line:
            break

    for line in lines[no+6:]:
        if 'TOTL' in line:
            _, value, sigma = float(info[0]), float(info[1]), float(info[3])
            total = value, sigma
            break

        info = line.split()
        group, value, sigma = float(info[0]), float(info[1]), float(info[2])
        flux.append(value)
        #flux[group] = value, sigma

    return total, flux

def getEnergyBoundary(out):

    energy = []
    lines = open(out).readlines()

    for no, line in enumerate(lines):
        if 'Energy Boundaries DEFAULT' in line:
            break
    for line in lines[no+4: ]:
        info = line.split()
        if len(info) == 1:
            break
        energy.append(float(info[1]))

    return energy


def collapseFlux(flux, groupStructure):
    '''
    '''
    collapsed = np.zeros(groupStructure)
    if groupStructure == 2:
        for i in range(len(flux)):
            if i < 214:
                collapsed[0] += flux[i]
            else:
                collapsed[1] += flux[i]
    return collapsed


def computeFourFactor(dir, PF=55, status='dry'):
    '''
    Compute the for factors using two group reaction rates 
    '''
    # quantities extracted from output files
    k, k_sigma, nubar, mfp, mfp_sigma = [], [], [], [], []
    total_absorption, total_fission, total_leakage = [], [], []

    # reaction rates read form .kmt files
    U235_absorption, U235_fission, U238_capture, U238_total = [], [], [], []

    # derived quantities
    eta, f, eps, p = [], [], [], []
    eta_sigma, f_sigma, eps_sigma, p_sigma = [], [], [], []

    thermal_leakage, fast_leakage = [], []
    fuel_thermal_absorption, fuel_fast_absorption = [], []

    pitchs = range(60, 92, 2)
    for no, pitch in enumerate(pitchs):
        case = 'VP55_v7.1-252_PF={}_pitch={}_en=9.6'.format(PF, pitch)
        if not os.path.isfile(dir + case + '.out'):
            case = 'VP55_v7.1-252_pebbles_PF={}_pitch={}_en=9.6'.format(PF, pitch)
        print('processing case {} .................'.format(case))

        # process output file
        outName = dir + case + '.out'
        keff, keff_sigma, runTime, mfp_, nubar_ ,\
        total_absorption_, total_fission_, total_leakage_ = kenoOut(outName)

        # process .kmt file
        reactionRates = getReaction(dir + case + '.kmt')

        # process reaction data from .out file
        outData = getKENOReactinData(outName)
        groupData = getKENOGroupReactinData(outName)

        k.append(keff)
        k_sigma.append(keff_sigma)

        nubar.append(nubar_[0])

        # total u235 absorption
        U235_absorption.append(reactionRates['92235'][27][-1][0])
        # thermal u235 fission
        U235_fission.append(reactionRates['92235'][18][1][0])
        # fast u238 capture
        U238_capture.append(reactionRates['92238'][102][0][0])
        # total u238 
        U238_total.append(reactionRates['92238'][1][-1][0])
        # total fission in system
        total_fission.append(total_fission_)
        # total leakage
        total_leakage.append(total_leakage_)
 
        # mean free path
        mfp.append(mfp_[0])
        mfp_sigma.append(mfp_[1])
       
        fuel_thermal_absorption.append(groupData['thermal']['fuel']['absorption'][0])
        fuel_fast_absorption.append(groupData['fast']['fuel']['absorption'][0])
        fast_leakage.append(groupData['fast']['total']['leakage'][0])
        thermal_leakage.append(groupData['thermal']['total']['leakage'][0])
        # thermal utilization factor
        f.append(
                 (groupData['thermal']['fuel']['absorption'][0]   +
                  groupData['thermal']['fuel']['fission'][0])     /
                 (groupData['thermal']['total']['absorption'][0]  + 
                  groupData['thermal']['total']['fission'][0])
                )
    
        f_sigma.append(np.sqrt( groupData['thermal']['fuel']['absorption'][1]**2
                              + groupData['thermal']['fuel']['fission'][1]**2 
                              + groupData['thermal']['total']['absorption'][1]**2
                              + groupData['thermal']['total']['fission'][1])**2
                      )

        # eta
        eta.append(nubar_[0]*(groupData['thermal']['fuel']['fission'][0]    +
                              groupData['fast']['fuel']['fission'][0])      /
                             (groupData['thermal']['fuel']['absorption'][0] + 
                              groupData['fast']['fuel']['absorption'][0])
                  )

        eta_sigma.append(np.sqrt(nubar_[1]**2 + (outData['1']['total_absorption'][1])**2
                         + 2*(total_fission_[1]**2)))

        # eps
        eps.append(
                   (groupData['fast']['fuel']['fission'][0] + groupData['thermal']['fuel']['fission'][0])/
                    groupData['thermal']['fuel']['fission'][0]
                   )

        eps_sigma.append(np.sqrt(reactionRates['92235'][18][-1][-1]**2 \
                         + reactionRates['92235'][18][1][-1]**2 + reactionRates['92238'][18][1][1]**2
                         +  reactionRates['92238'][18][-1][1]**2))

        # p = thermal_absorption/total_absorption
        # p.append(groupData['thermal']['total']['absorption'][0]/(groupData['thermal']['total']['absorption'][0] + 
        #                                                          groupData['fast']['total']['absorption'][0]))
        
        p_ = (groupData['thermal']['total']['absorption'][0] + 
              groupData['thermal']['total']['fission'][0]    +
              groupData['thermal']['total']['leakage'][0]) / \
             (groupData['fast']['total']['leakage'][0]       + 
              groupData['fast']['total']['absorption'][0]    +
              groupData['fast']['total']['fission'][0]       +
              groupData['thermal']['total']['fission'][0]    +
              groupData['thermal']['total']['leakage'][0]    +
              groupData['thermal']['total']['absorption'][0] )

        p.append(p_)

        sigma = np.sqrt(2*(groupData['thermal']['total']['absorption'][1]**2) +
                           groupData['fast']['total']['absorption'][1]**2)
        p_sigma.append(sigma)

    return {'k'  : (np.array(k), np.array(k_sigma)),
            'f'  : (np.array(f), np.array(f_sigma)),
            'p'  : (np.array(p), np.array(p_sigma)),
            'eps': (np.array(eps), np.array(eps_sigma)),
            'eta': (np.array(eta), np.array(eta_sigma)),
            'mfp': (np.array(mfp), np.array(mfp_sigma))
           }