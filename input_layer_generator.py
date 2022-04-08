# input: dimer_summary.dat (file with the dimer names, G_rel, E_exc)
# output: files MVR.out	MVS.out	QBB.out	QBF.out	QBS.out
# generate input layer for QBB model
infile = open('dimer_summary.dat','r')
conn = []
rox = [[0 for a in range(6)] for b in range(830)]
rbond = []
nl = 0
pi = -1
cis = str('c')
trans = str('t')
for line in infile:
    nl += 1
    pi += 1
    name = line.rstrip().split()[0]
    rconn = []
    if 'lin' in name:
        rconn.append(name[8:10])
        rox[pi][0] = name[13]
        rox[pi][1] = name[14]
        rox[pi][2] = name[15]
        rox[pi][3] = name[17]
        rox[pi][4] = name[18]
        rox[pi][5] = name[19]
        alfa = str(name[6])
        if cis in alfa:
           rbond.append('1')
        elif trans in alfa:
           rbond.append('2')
    elif 'cyc' in name:
        rbond.append('0')
        if name[6] > name[9]:
            rconn.append(name[9]+name[6])
        elif name[6] < name[9]:
            rconn.append(name[6]+name[9])
        elif name[6] == name[9]:
            rconn.append(name[6]+name[9])
        if name[7] > name[10]:
            rconn.append(name[10]+name[7])
        elif name[7] < name[10]:
            rconn.append(name[7]+name[10])
        elif name[7] == name[10]:
            rconn.append(name[7]+name[10])
        rox[pi][0] = name[14]
        rox[pi][1] = name[15]
        rox[pi][2] = name[16]
        rox[pi][3] = name[18]
        rox[pi][4] = name[19]
        rox[pi][5] = name[20]
    for k in range(len(rconn)):
        iflag = 0
        for j in range(len(conn)):
            if rconn[k] == conn[j]:
                iflag = 1
                break
        if iflag == 0:
            conn.append(rconn[k])
conn.sort()
for c in range(830):
    for d in range(6):
        if rox[c][d] == 'x':
           rox[c][d] = 0
molec = [[0 for x in range(len(conn))] for y in range(nl)]
infile.seek(0)
nl1 = -1
beta = 0
nom = [0 for x in range(nl)]
energia = [0 for z in range(nl)]
excitat = [0 for t in range(nl)]
for line in infile:
    nl1 += 1
    name = line.rstrip().split()[0]
    energy = line.rstrip().split()[1]
    excited = line.rstrip().split()[2]
    nom[beta] = name
    energia[beta] = energy
    excitat[beta] = excited
    rconn = []
    if 'lin' in name:
        rconn.append(name[8:10])
    elif 'cyc' in name:
        if name[6] > name[9]:
            rconn.append(name[9]+name[6])
        elif name[6] < name[9]:
            rconn.append(name[6]+name[9])
        elif name[6] == name[9]:
            rconn.append(name[6]+name[9])
        if name[7] > name[10]:
            rconn.append(name[10]+name[7])
        elif name[7] < name[10]:
            rconn.append(name[7]+name[10])
        elif name[7] == name[10]:
            rconn.append(name[7]+name[10])
    for k in range(len(rconn)):
        for j in range(len(conn)):
            if rconn[k] == conn[j]:
                if 'lin' in name:
                        molec[nl1][j] = 1
                elif 'cyc' in name:
                    if k == 0:
                        if name[6] > name[9]:
                            molec[nl1][j] = 2
                        elif name[6] < name[9]:
                            molec[nl1][j] = 1
                        elif name[6] == name[9]:
                            molec[nl1][j] = 1
                    if k == 1:
                        if name[7] > name[10]:
                            molec[nl1][j] = 2
                        elif name[7] < name[10]:
                            molec[nl1][j] = 1
                        elif name[7] == name[10]:
                            molec[nl1][j] = 1
                    if rconn[0] == rconn[1]:
                        molec[nl1][j] = 3
                        break
                        break
    beta += 1
infile.close()
outfile = open('QBB.out','w')
for k in range(nl+1):
    if k == 0:
        outfile.write('%s' % (''.join(str('Molec_name'))))
        outfile.write('%s' % (' '))
        for p in range(len(conn)):
            outfile.write('%s' % (''.join(str('Conn'))))
            outfile.write('%s' % ('_'))
            outfile.write('%s' % (''.join(str(conn[p]))))
            outfile.write('%s' % (' '))
        for t in range(6):
            position = t+1
            outfile.write('%s' % (''.join(str('Ox_position'))))
            outfile.write('%s' % ('_'))
            outfile.write('%s' % (''.join(str(position))))
            outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('Central_bond'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('G_rel'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('E_exc'))))
        outfile.write('%s\n' % (''))
    else:
        hola = k-1
        outfile.write('%s' % (''.join(str(nom[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(molec[hola][j]) for j in range(len(conn)))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(rox[hola][d]) for d in range(6))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(rbond[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(energia[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(excitat[hola]))))
        outfile.write('%s\n' % (''))
outfile.close()

# generate input layer for QBF model

infile = open('dimer_summary.dat','r')
conn = []
rox = [[0 for a in range(6)] for b in range(830)]
rbond = []
nl = 0
pi = -1
cis = str('c')
trans = str('t')
for line in infile:
    nl += 1
    pi += 1
    name = line.rstrip().split()[0]
    rconn = []
    if 'lin' in name:
        rconn.append(name[8])
        rconn.append(name[9])
        rox[pi][0] = name[13]
        rox[pi][1] = name[14]
        rox[pi][2] = name[15]
        rox[pi][3] = name[17]
        rox[pi][4] = name[18]
        rox[pi][5] = name[19]
        alfa = str(name[6])
        if cis in alfa:
           rbond.append('0')
        elif trans in alfa:
           rbond.append('1')
    elif 'cyc' in name:
        rconn.append(name[6]+name[7])
        if name[9] > name[10]:
            rbond.append('1')
            rconn.append(name[10]+name[9])
        elif name[9] < name[10]:
            rbond.append('0')
            rconn.append(name[9]+name[10])
        rox[pi][0] = name[14]
        rox[pi][1] = name[15]
        rox[pi][2] = name[16]
        rox[pi][3] = name[18]
        rox[pi][4] = name[19]
        rox[pi][5] = name[20]
    for k in range(len(rconn)):
        iflag = 0
        for j in range(len(conn)):
            if rconn[k] == conn[j]:
                iflag = 1
                break
        if iflag == 0:
            conn.append(rconn[k])
conn.sort()
for c in range(830):
    for d in range(6):
        if rox[c][d] == 'x':
           rox[c][d] = 0
molec = [[0 for x in range(len(conn))] for y in range(nl)]
infile.seek(0)
nl1 = -1
beta = 0
nom = [0 for x in range(nl)]
energia = [0 for z in range(nl)]
excitat = [0 for t in range(nl)]
for line in infile:
    nl1 += 1
    name = line.rstrip().split()[0]
    energy = line.rstrip().split()[1]
    excited = line.rstrip().split()[2]
    nom[beta] = name
    energia[beta] = energy
    excitat[beta] = excited
    rconn = []
    if 'lin' in name:
        rconn.append(name[8])
        rconn.append(name[9])
    elif 'cyc' in name:
        rconn.append(name[6]+name[7])
        if name[9] > name[10]:
            rconn.append(name[10]+name[9])
        elif name[9] < name[10]:
            rconn.append(name[9]+name[10])
    for k in range(len(rconn)):
        for j in range(len(conn)):
            if rconn[k] == conn[j]:
                molec[nl1][j] = 1
                break
                break
    beta += 1
infile.close()
outfile = open('QBF.out','w')
for k in range(nl+1):
    if k == 0:
        outfile.write('%s' % (''.join(str('Molec_name'))))
        outfile.write('%s' % (' '))
        for p in range(len(conn)):
            outfile.write('%s' % (''.join(str('Conn'))))
            outfile.write('%s' % ('_'))
            outfile.write('%s' % (''.join(str(conn[p]))))
            outfile.write('%s' % (' '))
        for t in range(6):
            position = t+1
            outfile.write('%s' % (''.join(str('Ox_position'))))
            outfile.write('%s' % ('_'))
            outfile.write('%s' % (''.join(str(position))))
            outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('Central_bond'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('G_rel'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('E_exc'))))
        outfile.write('%s\n' % (''))
    else:
        hola = k-1
        outfile.write('%s' % (''.join(str(nom[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(molec[hola][j]) for j in range(len(conn)))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(rox[hola][d]) for d in range(6))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(rbond[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(energia[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(excitat[hola]))))
        outfile.write('%s\n' % (''))
outfile.close()

# generate input layer for QBS model

infile = open('dimer_summary.dat','r')
conn = []
conn2 = []
rox = [[0 for a in range(6)] for b in range(830)]
rbond = []
nl = 0
pi = -1
cis = str('c')
trans = str('t')
for line in infile:
    nl += 1
    pi += 1
    name = line.rstrip().split()[0]
    rconn = []
    rconn2 = []
    if 'lin' in name:
        rconn.append(name[8])
        rconn2.append(name[9])
        rox[pi][0] = name[13]
        rox[pi][1] = name[14]
        rox[pi][2] = name[15]
        rox[pi][3] = name[17]
        rox[pi][4] = name[18]
        rox[pi][5] = name[19]
        alfa = str(name[6])
        if cis in alfa:
           rbond.append('0')
        elif trans in alfa:
           rbond.append('1')
    elif 'cyc' in name:
        rconn.append(name[6])
        rconn.append(name[7])
        rconn2.append(name[9])
        rconn2.append(name[10])
        if name[9] > name[10]:
            rbond.append('1')
        elif name[9] < name[10]:
            rbond.append('0')
        rox[pi][0] = name[14]
        rox[pi][1] = name[15]
        rox[pi][2] = name[16]
        rox[pi][3] = name[18]
        rox[pi][4] = name[19]
        rox[pi][5] = name[20]
    for k in range(len(rconn)):
        iflag = 0
        for j in range(len(conn)):
            if rconn[k] == conn[j]:
                iflag = 1
                break
        if iflag == 0:
            conn.append(rconn[k])
    for k2 in range(len(rconn2)):
        iflag = 0
        for j2 in range(len(conn2)):
            if rconn2[k2] == conn2[j2]:
                iflag = 1
                break
        if iflag == 0:
            conn2.append(rconn2[k2])
conn.sort()
conn2.sort()
for c in range(830):
    for d in range(6):
        if rox[c][d] == 'x':
           rox[c][d] = 0
molec = [[0 for x in range(len(conn))] for y in range(nl)]
molec2 = [[0 for x in range(len(conn2))] for y in range(nl)]
infile.seek(0)
nl1 = -1
beta = 0
nom = [0 for x in range(nl)]
energia = [0 for z in range(nl)]
excitat = [0 for t in range(nl)]
for line in infile:
    nl1 += 1
    name = line.rstrip().split()[0]
    energy = line.rstrip().split()[1]
    excited = line.rstrip().split()[2]
    nom[beta] = name
    energia[beta] = energy
    excitat[beta] = excited
    rconn = []
    rconn2 = []
    if 'lin' in name:
        rconn.append(name[8])
        rconn2.append(name[9])
    elif 'cyc' in name:
        rconn.append(name[6])
        rconn.append(name[7])
        rconn2.append(name[9])
        rconn2.append(name[10])
    for k in range(len(rconn)):
        for j in range(len(conn)):
            if rconn[k] == conn[j]:
                if 'lin' in name:
                    molec[nl1][j] = 1
                if 'cyc' in name:
                    molec[nl1][j] = 2
                break
                break
    for k2 in range(len(rconn2)):
        for j2 in range(len(conn2)):
            if rconn2[k2] == conn2[j2]:
                if 'lin' in name:
                    molec2[nl1][j2] = 1
                if 'cyc' in name:
                    molec2[nl1][j2] = 2
                break
                break
    beta += 1
infile.close()
outfile = open('QBS.out','w')
for k in range(nl+1):
    if k == 0:
        outfile.write('%s' % (''.join(str('Molec_name'))))
        outfile.write('%s' % (' '))
        for p in range(len(conn)):
            outfile.write('%s' % (''.join(str('Conn'))))
            outfile.write('%s' % ('_'))
            outfile.write('%s' % (''.join(str(conn[p]))))
            outfile.write('%s' % (' '))
        for p2 in range(len(conn2)):
            outfile.write('%s' % (''.join(str('Conn2'))))
            outfile.write('%s' % ('_'))
            outfile.write('%s' % (''.join(str(conn2[p2]))))
            outfile.write('%s' % (' '))
        for t in range(6):
            position = t+1
            outfile.write('%s' % (''.join(str('Ox_position'))))
            outfile.write('%s' % ('_'))
            outfile.write('%s' % (''.join(str(position))))
            outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('Central_bond'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('G_rel'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('E_exc'))))
        outfile.write('%s\n' % (''))
    else:
        hola = k-1
        outfile.write('%s' % (''.join(str(nom[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(molec[hola][j]) for j in range(len(conn)))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(molec2[hola][j2]) for j2 in range(len(conn2)))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(rox[hola][d]) for d in range(6))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(rbond[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(energia[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(excitat[hola]))))
        outfile.write('%s\n' % (''))
outfile.close()

# generate input layer for MVR model

infile = open('dimer_summary.dat','r')
conn = []
rox = [[0 for a in range(6)] for b in range(830)]
rbond = []
nl = 0
pi = -1
cis = str('c')
trans = str('t')
ox = []
tox = []
for line in infile:
    nl += 1
    pi += 1
    name = line.rstrip().split()[0]
    rconn = []
    if 'lin' in name:
        rconn.append(name[8:10])
        rox[pi][0] = name[13]
        rox[pi][1] = name[14]
        rox[pi][2] = name[15]
        rox[pi][3] = name[17]
        rox[pi][4] = name[18]
        rox[pi][5] = name[19]
        alfa = str(name[6])
        if cis in alfa:
           rbond.append('1')
        elif trans in alfa:
           rbond.append('2')
    elif 'cyc' in name:
        rconn.append(name[6:8]+name[9:11])
        rox[pi][0] = name[14]
        rox[pi][1] = name[15]
        rox[pi][2] = name[16]
        rox[pi][3] = name[18]
        rox[pi][4] = name[19]
        rox[pi][5] = name[20]
        rbond.append('0')
    for k in range(len(rconn)):
        iflag = 0
        for j in range(len(conn)):
            if rconn[k] == conn[j]:
                iflag = 1
                break
        if iflag == 0:
            conn.append(rconn[k])
    for d in range(6):
        if rox[pi][d] == 'x':
           rox[pi][d] = '0'
    tox.append(rox[pi][0]+rox[pi][1]+rox[pi][2]+rox[pi][3]+rox[pi][4]+rox[pi][5])
    for e in range(len(tox)):
        fflag = 0
        for z in range(len(ox)):
            if tox[e] == ox[z]:
                fflag = 1
                break
        if fflag == 0:
            ox.append(tox[e])
conn.sort()
ox.sort()
molec = [0 for y in range(nl)]
gox = [0 for v in range(nl)]
infile.seek(0)
nl1 = -1
beta = 0
nom = [0 for x in range(nl)]
energia = [0 for z in range(nl)]
excitat = [0 for t in range(nl)]
boom = -1
yey = [0 for g in range(len(conn))]
k = 0
j = 0
yea = [0 for s in range(len(ox))]
for g in range(len(conn)):
    yey[g] = boom+1
    boom = boom+1
boom = -1
for s in range(len(ox)):
    yea[s] = boom+1
    boom = boom+1
for line in infile:
    nl1 += 1
    name = line.rstrip().split()[0]
    energy = line.rstrip().split()[1]
    excited = line.rstrip().split()[2]
    nom[beta] = name
    energia[beta] = energy
    excitat[beta] = excited
    rconn = []
    tox = []
    if 'lin' in name:
        rconn.append(name[8:10])
        tox.append(name[13].replace('x','0')+name[14].replace('x','0')+name[15].replace('x','0')+name[17].replace('x','0')+name[18].replace('x','0')+name[19].replace('x','0'))
    elif 'cyc' in name:
        rconn.append(name[6:8]+name[9:11])
        tox.append(name[14].replace('x','0')+name[15].replace('x','0')+name[16].replace('x','0')+name[18].replace('x','0')+name[19].replace('x','0')+name[20].replace('x','0'))
    for k in range(len(rconn)):
        for j in range(len(conn)):
            if rconn[k] == conn[j]:
                molec[nl1] = yey[j]
                break
    for e in range(len(tox)):
        for z in range(len(ox)):
            if tox[e] == ox[z]:
                gox[nl1] = yea[z]
    beta += 1
infile.close()
outfile = open('MVR.out','w')
for k in range(nl+1):
    if k == 0:
        outfile.write('%s' % (''.join(str('Molec_name'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('Conn'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('Ox_position'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('Double_bond'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('G_rel'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('E_exc'))))
        outfile.write('%s\n' % (''))
    else:
        hola = k-1
        outfile.write('%s' % (''.join(str(nom[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(molec[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(gox[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(rbond[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(energia[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(excitat[hola]))))
        outfile.write('%s\n' % (''))
outfile.close()

# generate input layer for MVS model

infile = open('dimer_summary.dat','r')
conn = []
rox = [[0 for a in range(6)] for b in range(830)]
rbond = []
nl = 0
pi = -1
cis = str('c')
trans = str('t')
for line in infile:
    nl += 1
    pi += 1
    name = line.rstrip().split()[0]
    rconn = []
    if 'lin' in name:
        rconn.append(name[8])
        rconn.append(name[9])
        rox[pi][0] = name[13]
        rox[pi][1] = name[14]
        rox[pi][2] = name[15]
        rox[pi][3] = name[17]
        rox[pi][4] = name[18]
        rox[pi][5] = name[19]
        alfa = str(name[6])
        if cis in alfa:
           rbond.append('0')
        elif trans in alfa:
           rbond.append('1')
    elif 'cyc' in name:
        rconn.append(name[6])
        rconn.append(name[7])
        rconn.append(name[9])
        rconn.append(name[10])
        if name[9] > name[10]:
            rbond.append('1')
        elif name[9] < name[10]:
            rbond.append('0')
        rox[pi][0] = name[14]
        rox[pi][1] = name[15]
        rox[pi][2] = name[16]
        rox[pi][3] = name[18]
        rox[pi][4] = name[19]
        rox[pi][5] = name[20]
    for k in range(len(rconn)):
        iflag = 0
        for j in range(len(conn)):
            if rconn[k] == conn[j]:
                iflag = 1
                break
        if iflag == 0:
            conn.append(rconn[k])
conn.sort()
for c in range(830):
    for d in range(6):
        if rox[c][d] == 'x':
           rox[c][d] = 0
roxfinal = [[0 for a in range(3)] for b in range(830)]
molec = [[0 for x in range(len(conn))] for y in range(nl)]
infile.seek(0)
nl1 = -1
beta = 0
nom = [0 for x in range(nl)]
energia = [0 for z in range(nl)]
excitat = [0 for t in range(nl)]
for line in infile:
    nl1 += 1
    name = line.rstrip().split()[0]
    energy = line.rstrip().split()[1]
    excited = line.rstrip().split()[2]
    nom[beta] = name
    energia[beta] = energy
    excitat[beta] = excited
    rconn = []
    if 'lin' in name:
        rconn.append(name[8])
        rconn.append(name[9])
    elif 'cyc' in name:
        rconn.append(name[6])
        rconn.append(name[7])
        rconn.append(name[9])
        rconn.append(name[10])
    for k in range(len(rconn)):
        for j in range(len(conn)):
            if 'lin' in name:
                if rconn[k] == conn[j]:
                    if rconn[0] == rconn[1]:
                        molec[nl1][j] = 3
                    elif k==0:
                        molec[nl1][j] = 1
                    elif k==1:
                        molec[nl1][j] = 2
            if 'cyc' in name:
                if rconn[k] == conn[j]:
                    if k==0 or k==1:
                        if rconn[k] == rconn[2] or rconn[k] == rconn[3]:
                            molec[nl1][j] = 3
                        else:
                            molec[nl1][j] = 1
                        break
                    elif k==2 or k==3:
                        if rconn[k] == rconn[0] or rconn[k] == rconn[1]:
                            molec[nl1][j] = 3
                        else:
                            molec[nl1][j] = 2
                        break
    for f in range(6):
        if rox[nl1][f] == '1':
            if f==0:
                if rox[nl1][f] == rox[nl1][3]:
                    roxfinal[nl1][f] = 3
                else:
                    roxfinal[nl1][f] = 1
            if f==1:
                if rox[nl1][f] == rox[nl1][4]:
                    roxfinal[nl1][f] = 3
                else:
                    roxfinal[nl1][f] = 1
            if f==2:
                if rox[nl1][f] == rox[nl1][5]:
                    roxfinal[nl1][f] = 3
                else:
                    roxfinal[nl1][f] = 1
            if f==3:
                if rox[nl1][f] == rox[nl1][0]:
                    roxfinal[nl1][f-3] = 3
                else:
                    roxfinal[nl1][f-3] = 2
            if f==4:
                if rox[nl1][f] == rox[nl1][1]:
                    roxfinal[nl1][f-3] = 3
                else:
                    roxfinal[nl1][f-3] = 2
            if f==5:
                if rox[nl1][f] == rox[nl1][2]:
                    roxfinal[nl1][f-3] = 3
                else:
                    roxfinal[nl1][f-3] = 2
    beta += 1
infile.close()
outfile = open('MVS.out','w')
for k in range(nl+1):
    if k == 0:
        outfile.write('%s' % (''.join(str('Molec_name'))))
        outfile.write('%s' % (' '))
        for p in range(len(conn)):
            outfile.write('%s' % (''.join(str('Conn'))))
            outfile.write('%s' % ('_'))
            outfile.write('%s' % (''.join(str(conn[p]))))
            outfile.write('%s' % (' '))
        for t in range(3):
            position = t+1
            outfile.write('%s' % (''.join(str('Ox_position'))))
            outfile.write('%s' % ('_'))
            outfile.write('%s' % (''.join(str(position))))
            outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('Central_bond'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('G_rel'))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str('E_exc'))))
        outfile.write('%s\n' % (''))
    else:
        hola = k-1
        outfile.write('%s' % (''.join(str(nom[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(molec[hola][j]) for j in range(len(conn)))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(roxfinal[hola][d]) for d in range(3))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (' '.join(str(rbond[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(energia[hola]))))
        outfile.write('%s' % (' '))
        outfile.write('%s' % (''.join(str(excitat[hola]))))
        outfile.write('%s\n' % (''))
outfile.close()

