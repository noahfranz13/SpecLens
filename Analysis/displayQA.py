'''
Simple script to write an html page to display 
quality analysis (QA) pages from fastspecfit
'''
import os, sys, glob
import argparse

def main():

    p = argparse.ArgumentParser()
    p.add_argument('--dir', help='Directory to create html page in', required=True)
    args = p.parse_args()

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
    <th>QA Plot</th>
    </tr>
    '''
    
    with open(idxfile, 'w') as f:
        f.write(initStr)
        for png in pngs:
            f.write('<tr>\n')
            f.write(f'<td><img src={os.path.split(png)[1]} height="60%" width="75%"></td>\n')
            f.write(f'</tr>')
        f.write('</table>')
        f.write('</html></body>')
            
if __name__ == '__main__':
    sys.exit(main())
