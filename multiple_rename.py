import os
# Function to rename multiple files
def main():
   i = 1
   path="./"
   for filename1 in os.listdir(path):
      if 'Spectra_Condor.zip.0' in filename1:
        filename2 = filename1.replace( '.zip', '.mat.zip' )
        print( i, filename1, '->', filename2 )
        # rename() function will
        # rename all the files
        os.rename( filename1, filename2 )
        i += 1
# Driver Code
if __name__ == '__main__':
   # Calling main() function
   main()