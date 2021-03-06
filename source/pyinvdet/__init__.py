from pkg_resources import get_distribution, DistributionNotFound
import os.path

try:
    _dist = get_distribution("invdet")
    dist_loc = os.path.realpath( os.path.normcase(_dist.location) )
    here = os.path.realpath( os.path.normcase(__file__) )
    if not here.startswith(os.path.join(dist_loc, 'invdet')):
        raise DistributionNotFound
except DistributionNotFound:
    __version__ = 'Please install this project with pip'
else:
    __version__ = _dist.version
