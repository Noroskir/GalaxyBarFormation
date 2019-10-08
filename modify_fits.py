import numpy as np
from astropy.io import fits

hdul = fits.open("data/test.fits")

data = hdul[1].data


def remove_row(data, index):
    for i in range(len(data)):
        data[i] = np.delete(data[i], index)
    return data


def investigate(data, index):
    for i in index:
        print("Index: {}".format(i))
        print("Bar: {}\t Disk: {}\t Bulge: {}".format(
            data[0][i], data[1][i], data[2][i]))


def get_type(name):
    dat = fits.open('data/sample/CALIFA_2_MS_class.fits')[1].data
    ind = np.where(dat['REALNAME'] == name)[0]
    if len(ind) == 0:
        dat = fits.open('data/sample/CALIFA_2_ES_class.fits')[1].data
        ind = np.where(dat['realname'] == name)[0]
    return (dat['hubtyp'][ind] + dat['hubsubtyp'][ind])[0]


def clear_data(data, g, name):
    ind = np.where(g == name)[0]
    typ = get_type(name)
    if typ[0] == 'E':
        if not disk[ind[0]]:
            ind = ind[0]
        else:
            ind = ind[1]
    elif typ[0] == 'S' and data[0][ind[0]] == data[0][ind[1]]:
        if disk[ind[0]]:
            ind = ind[0]
        else:
            ind = ind[1]
    elif typ[0] == 'S' and data[0][ind[0]] != data[0][ind[1]]:
        if not data[0][ind[0]]:
            ind = ind[0]
        else:
            ind = ind[1]
    else:
        print('In clear_data case not implemented!')
        print('Investigating:')
        investigate(data, ind)
    g = np.delete(g, ind)
    data = remove_row(data, ind)
    return g, data


g = data['GALAXY'][0]
bar = data['BAR'][0]
disk = data['DISK'][0]
bulge = data['BULGE'][0]
dbreak = data['DBREAK'][0]
redshift = data['REDSHIFT'][0]
mue = data['MUE_R'][0]
re = data['RE_R'][0]
n = data['N_R'][0]
mu0b = data['MU0B_R'][0]
nbar = data['NBAR_R'][0]
rbar = data['RBAR_R'][0]
mu0 = data['MU0_R'][0]
hi = data['HI_R'][0]
ho = data['HO_R'][0]
rbreak = data['RBREAK_R'][0]
babar = data['BABAR_R'][0]

ndata = ['BAR', 'DISK', 'BULGE', 'DBREAK', 'REDSHIFT', 'MUE_R',
         'RE_R', 'N_R', 'MU0B_R', 'NBAR_R', 'RBAR_R', 'MU0_R',
         'HI_R', 'HO_R', 'RBREAK_R', 'BABAR_R']
form = ['I', 'I', 'I', 'I', 'E', 'D', 'D', 'D', 'D', 'D',
        'D', 'D', 'D', 'D', 'D', 'D']
data = [bar, disk, bulge, dbreak, redshift, mue, re, n, mu0b, nbar,
        rbar, mu0, hi, ho, rbreak, babar]

# bulge only or with disk, classified as bulge
ind = np.where(g == 'UGC00029')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# classified as elliptical, remove disk
ind = np.where(g == 'NGC0499')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# classified as S0, remove bulk only
ind = np.where(g == 'NGC0932')[0]
if disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# elliptical, remove disk
ind = np.where(g == 'NGC0962')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# classified as S0, remove bulge only
ind = np.where(g == 'UGC02099')[0]
if disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# elliptical, remove disk
ind = np.where(g == 'NGC1060')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# S0 remove bulge only
ind = np.where(g == 'NGC1167')[0]
if disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# E2
ind = np.where(g == 'NGC2513')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# E
ind = np.where(g == 'NGC3615')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# E
ind = np.where(g == 'NGC4816')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# E
ind = np.where(g == 'NGC6021')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# S
ind = np.where(g == 'NGC6173')[0]
if disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# E
ind = np.where(g == 'NGC7194')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# E
ind = np.where(g == 'NGC7562')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# S0
ind = np.where(g == 'NGC7683')[0]
if disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# E
ind = np.where(g == 'UGC10693')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# S0
ind = np.where(g == 'UGC10905')[0]
if disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)

# E
ind = np.where(g == 'UGC12127')[0]
if not disk[ind[0]]:
    ind = ind[0]
else:
    ind = ind[1]
g = np.delete(g, ind)
data = remove_row(data, ind)


g, data = clear_data(data, g, 'NGC1132')
g, data = clear_data(data, g, 'NGC5423')
g, data = clear_data(data, g, 'NGC5580')
g, data = clear_data(data, g, 'NGC5687')
g, data = clear_data(data, g, 'IC3065')
g, data = clear_data(data, g, 'IC3652')
g, data = clear_data(data, g, 'NGC6023')
g, data = clear_data(data, g, 'IC1602')
g, data = clear_data(data, g, 'UGC03960')
g, data = clear_data(data, g, 'MCG+07-17-002')
g, data = clear_data(data, g, 'NGC2484')
g, data = clear_data(data, g, 'IC2378')
#g, data = clear_data(data, g, 'UGC04458')
g, data = clear_data(data, g, 'IC2402')
#g, data = clear_data(data, g, '2MASXJ09065870')
g, data = clear_data(data, g, 'CGCG429-012')
g, data = clear_data(data, g, 'NGC0426')
g, data = clear_data(data, g, 'NGC0472')


# create fits file
hdr = hdul[0].header
hdu1 = fits.PrimaryHDU(header=hdr)


for i in range(len(data)):
    data[i] = fits.Column(name=ndata[i], format=form[i], array=data[i])

g = fits.Column(name='GALAXY', format='6600A', array=g)

data = [g] + data

cols = fits.ColDefs(data)
hdu2 = fits.BinTableHDU.from_columns(cols)

hdul = fits.HDUList([hdu1, hdu2])
hdul.writeto('data/cut_photometrics.fits')
