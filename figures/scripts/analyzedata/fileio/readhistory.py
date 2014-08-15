import cPickle;
import os;

class histReader:
    floaded = False;
    histpaths = [];     # Previous paths
    f = '';
    maxhist = -1;
    
    def __init__(self, filename, max_hist):
        self.fname = filename;
        
        if max_hist >= 1:
            self.maxhist = max_hist;  # Update the max_hist.
    
        if os.path.exists(filename):
            try:
                self.reload();
            except:
                raise;
        
    def reload(self):
        if os.path.exists(self.fname):
            try:
                self.f = open(self.fname, 'r');
            
                self.floaded = True;
            
                self.maxhist = cPickle.load(self.f);
                self.histpaths = cPickle.load(self.f);
            
                self.f.close();
            except:
                raise;
        else:
            raise Exception('File does not exist.');
    
    def update(self):
        # Updates the file.
        if(self.maxhist < 1):
            self.maxhist = 20;        
        
        # Truncate the history path.
        histpaths = self.histpaths;
        del histpaths[self.maxhist:];
            
        try:           
            self.f = open(self.fname, 'w+');
            
            cPickle.dump(self.maxhist, self.f);
            cPickle.dump(histpaths, self.f);
            
            self.f.close();
        except:
            raise;        
    
    def addpath(self, path):
        # Updates the history path
        if path == '':
            raise Exception('Empty string passed to histpath.');
        
        try:
            self.histpaths.remove(path) # No problem if it's an error.
        except:
            pass;
        
        self.histpaths.insert(0, path);
        
        
        