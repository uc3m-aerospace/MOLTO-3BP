def queryFunc():

    print('Mission sequence simulation')

    # Iteration variables preallocated
    seq = {'it': 1}

    while 1:

        print('\nParameters for Orbit ' + str(seq['it']))

        seq['type' + str(seq['it'])] = input('Orbit type [HL, LY]: ')

        if seq['type' + str(seq['it'])] == 'LY':
            Ax = input('Amplitude (Ax, nondimensional): ')
            LP = input('Lagrange Point [1, 2]: ')

            seq['Ax' + str(seq['it'])] = float(Ax)
            seq['LP' + str(seq['it'])] = int(LP)

        elif seq['type' + str(seq['it'])] == 'HL':
            Az  = input('Amplitude (Az, dimensional [km]): ')
            phi = input('Azimuthal rotation [rads]: ')
            m   = input('North/South variants of the orbit [1, 3]: ')
            LP  = input('Lagrange Point [1, 2]: ')

            seq['Az' + str(seq['it'])]  = float(Az)
            seq['phi' + str(seq['it'])] = float(phi)
            seq['m' + str(seq['it'])]   = int(m)
            seq['LP' + str(seq['it'])]  = int(LP)

        else:
            raise Exception('Manifolds_Seq:typeError.'+\
                '    The type selected is not valid [\'LY\'][\'HL\']!')

        key = input('\nWould you like to add a section? (1 == yes, 0 == no): ')

        if not int(key):
            break

        print('\nParameters for Section ' + str(seq['it']))
        ang = input('Angle between section and +X semiplane [degrees]: ')

        seq['ang' + str(seq['it'])] = int(ang)

        key = input('\nWould you like to add another orbit? (1 == yes, 0 == no): ')

        if not int(key):
            break

        seq['it'] = seq['it'] +1

    return seq

def queryFunctext(Input):

    print('Mission sequence simulation')

    # Iteration variables preallocated
    seq  = {'it': 1}
    fid  = open(Input['file'], 'r')
    data = fid.readlines()
    refdata = []
    for i in range(len(data)):
        if data[i].split()[0][0] != '#':
            refdata.append(i)

    data = [data[j] for j in refdata]
    line = 0

    while 1:

        seq['type' + str(seq['it'])] = data[line].split()[-1]

        if seq['type' + str(seq['it'])] == 'LY':
            Ax = data[line +1].split()[-1]
            LP = data[line +2].split()[-1]

            seq['Ax' + str(seq['it'])] = float(Ax)
            seq['LP' + str(seq['it'])] = int(LP)

            if len(data) == line +3:
                break

            ang = data[line +3].split()[-1]

            seq['ang' + str(seq['it'])] = int(ang)

            if len(data) == line +4:
                break

            line = line +4

        elif seq['type' + str(seq['it'])] == 'HL':
            Az  = data[line +1].split()[-1]
            phi = data[line +2].split()[-1]
            m   = data[line +3].split()[-1]
            LP  = data[line +4].split()[-1]

            seq['Az' + str(seq['it'])]  = float(Az)
            seq['phi' + str(seq['it'])] = float(phi)
            seq['m' + str(seq['it'])]   = int(m)
            seq['LP' + str(seq['it'])]  = int(LP)

            if len(data) == line +5:
                break

            ang = data[line +5].split()[-1]

            seq['ang' + str(seq['it'])] = int(ang)

            if len(data) == line +6:
                break

            line = line +6

        else:
            raise Exception('Manifolds_Seq:typeError.'+\
                '    The type selected is not valid [\'LY\'][\'HL\']!')

        seq['it'] = seq['it'] +1

    fid.close()

    return seq
