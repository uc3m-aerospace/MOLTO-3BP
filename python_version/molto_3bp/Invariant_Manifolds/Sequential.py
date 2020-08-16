def queryFunc():

    print('Mission sequence simulation')

    # Iteration variables preallocated
    seq = {'it': 1}

    while 1:

        print('\nParameters for Orbit ' + str(seq['it']))
        Ax = input('Amplitude (Ax, nondimensional): ')
        LP = input('Lagrange Point [1, 2]: ')

        seq['Ax' + str(seq['it'])] = float(Ax)
        seq['LP' + str(seq['it'])] = int(LP)

        key = input('\nWould you like to add a section? (1 == yes, 0 == no): ')

        if not int(key):
            break

        print('\nParameters for Section ' + str(seq['it']))
        ang = input('Angle between section and +X semiplane: ')

        seq['ang' + str(seq['it'])] = int(ang)

        key = input('\nWould you like to add another orbit? (1 == yes, 0 == no): ')

        if not int(key):
            break

        seq['it'] = seq['it'] +1

    return seq
