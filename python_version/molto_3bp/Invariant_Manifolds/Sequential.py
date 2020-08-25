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

    while 1:

        Ax = data[3*seq['it']-3].split()[-1]
        LP = data[3*seq['it']-2].split()[-1]

        seq['Ax' + str(seq['it'])] = float(Ax)
        seq['LP' + str(seq['it'])] = int(LP)

        if len(data) == 3*seq['it']-1:
            break

        ang = data[3*seq['it']-1].split()[-1]

        seq['ang' + str(seq['it'])] = int(ang)

        if len(data) == 3*seq['it']:
            break

        seq['it'] = seq['it'] +1

    fid.close()

    return seq
