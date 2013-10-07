from os import listdir
from os.path import isfile, join

location = '/scratch/garbetsp/varsens/samples'

onlyfiles = [f for f in listdir(location) if isfile(join(location,f)) ]

objs  = [i for i in onlyfiles if i[0:9] == 'objective']

jobs = [int(f.split('-')[1].split('.')[0]) for f in objs]

m = max(jobs)


print("Range is 1 to ", m)

full_set = set(xrange(1, m))

print("Missing: ", sorted(list(full_set - set(jobs))))

