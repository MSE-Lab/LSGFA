from modules.pfam import *

hits = Hits('/var/folders/nh/9t3m8tz15cl8z_jd273ryx800000gn/T/out_yv6lmzgr')
for i in hits.ana_relations():
    print(i)
