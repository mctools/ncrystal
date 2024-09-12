#!/usr/bin/env python3

#import time
#time.sleep(3)
def sec(title):
    print()
    print ('='*100)
    print ('='*100)
    print(f'====  {title.center(88)}  ====')
    print ('='*100)
    print ('='*100)
    print()

def main():
    data = """NCMAT v1
@CELL
    lengths 4.04958 4.04958 4.04958
    angles 90. 90. 90.
@ATOMPOSITIONS
    Al 0. 0.5 0.5
    Al 0. 0. 0.
    Al 0.5 0.5 0.
    Al 0.5 0. 0.5
@DEBYETEMPERATURE
    Al   410.3542
"""

    import NCrystal as NC
    NC.removeAllDataSources()
    NC.registerInMemoryFileData('Al_nosg.ncmat',data)
    NC.registerInMemoryFileData('Al.ncmat',data+'@SPACEGROUP\n    225\n')
    for verbose in (None,0,1,2,3,'sc'):
        for fn in ('Al_nosg.ncmat','Al.ncmat'):
            if verbose is None:
                sec(f'Contents of "{fn}"')
                print(NC.createTextData(fn).rawData)
                continue
            cfg=f'{fn};dcutoff=0.4'
            if verbose=='sc':
                sec(f'Scatter-dump of "{cfg}"')
                NC.createScatter(cfg).dump()
            else:
                sec(f'Info-dump of "{cfg}" with verbosity lvl {verbose}')
                NC.createInfo(cfg).dump(verbose=verbose)

if __name__=='__main__':
    main()
