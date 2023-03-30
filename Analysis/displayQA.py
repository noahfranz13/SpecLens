'''
Simple script to write an html page to display 
quality analysis (QA) pages from fastspecfit
'''
import os, sys, glob
import argparse
import fitsio
from astropy.table import Table

def main():

    p = argparse.ArgumentParser()
    p.add_argument('--dir', help='Directory to create html page in', required=True)
    args = p.parse_args()

    # read in astropy tables to populate table with more info
    ironMatchesPath = os.path.join(os.path.split(args.dir)[0], 'iron-matches.fits')
    ironMatches = Table(fitsio.read(ironMatchesPath))

    # read in source redrock fit to only add plots with a good fit to the html page
    rrSourceFile = os.path.join(os.path.split(args.dir)[0], 'redrock-source.fits')
    rrSource = Table(fitsio.read(rrSourceFile))
    
    pngs = glob.glob(os.path.join(args.dir, '*.png'))
    
    idxfile = os.path.join(args.dir, 'index.html')

    initStr = ''' 
    <html><body>
    <h1>SpecLens fastspecfit QA Plots</h1>
    <style type="text/css">
    table, td, th {padding: 5px; text-align: center; border: 1px solid black;}
    p {display: inline-block;;}
    </style>
    <table>
    <tr>
    <th>Target Info</th>
    <th>QA Plot</th>
    </tr>
    '''
    
    with open(idxfile, 'w') as f:
        f.write(initStr)
        ii = 0
        for png in pngs:
            tid = png.split('.')[0].split('-')[-1]

            # skip files that have a bad fit
            I = rrSource['TARGETID'] == int(tid)
            rr = rrSource[I]
            if rr['ZWARN'] != 0:
                continue
            ii += 1
            flag = ironMatches['TARGETID'] == int(tid)
            iflag = ironMatches[flag]

            ra = iflag['TARGET_RA'][0]
            dec = iflag['TARGET_DEC'][0]
            
            col0 = f'''
            <td>
            Target ID: {tid}\n
            <br>RA: {ra}\n
            <br>Dec: {dec}\n
            <br><a href="http://legacysurvey.org/viewer-dev/?ra={ra}&dec={dec}&zoom=14&layer=ls-dr9&sga&sga-parent" target="_blank">Legacy Link</a>
            </td>\n
            '''      
            
            f.write('<tr>\n')
            f.write(col0)
            f.write(f'<td><img src={os.path.split(png)[1]} height="60%" width="75%"></td>\n')
            f.write(f'</tr>\n')
        f.write('</table>')
        f.write('</html></body>')

    print(f'QA Plots Displayed for {ii} Targets, the rest have ZWARN > 0!')
        
if __name__ == '__main__':
    sys.exit(main())
