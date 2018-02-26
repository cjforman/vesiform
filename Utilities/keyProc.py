'''
Created on 14 Dec 2017

@author: chris
'''
import sys
from fileIO import readTextFile

class keyProc():
    ''' Class constructs a dictionary called params which stores values from a text based parameter list consisting of lines in the format of:
    
    type <name> <value>
    
    Type can be float, string or int. <name> is the key and <value> is the value in the dictionary. E.G
    
    float Angle 23.5
    
    is equivalent to: params['Angle'] = 23.5 
    
    It is posssible to append values to a list of the same type using 
    
    stringlist, intlist, floatlist, stringlistbuild
    
    '''
    
    def __init__(self, paramFilename):
        
        self.paramFilename = paramFilename
        
        # set up the parameter dictionary if it hasn't been done already 
        self.params = {}

        # set up the error reporting variable
        self.noLoadErrors = True 
    
        # load the parameters
        self.loadParams()
        
        # report on initialisation process
        if self.noLoadErrors==True:
            print "Parameters loaded and initialised successfully"
        else:
            print "Error loading parameters"
            self.dumpParams()
            sys.exit()
    
    def loadParams(self):
        # populate the parameter dictionary from the input data function
        self.GetInputData()
        
        if self.noLoadErrors:
            self.parseParams()
        
        if self.noLoadErrors:
            self.initialiseParameters()
    
           
    def initialiseParameters(self):
        # Overload to perform preprocessing based on loaded parameters
        self.noLoadErrors = True
  
    def GetInputData(self):
        self.lines = readTextFile(self.paramFilename)
        if self.lines==False:
            self.noLoadErrors = False
    
    def getParam(self, paramString):
        # function to load a param and exit gracefully with an error report if it fails.
        # in this way a list of unspecified but required parameters is generated.
        paramVal = 'NULL'
        try:
            paramVal = self.params[paramString]
        except KeyError:
            print("Key error: Unknown key", paramString)
            self.noLoadErrors = False
            pass
        return paramVal
        
    def parseParams(self):
        noErrors = True
        for param in self.lines:
            if ',' in param:
                noErrors=False
                print("Use spaces to separate values not commas. Thanks.")
                break

            if noErrors:
                tokens = param.split()
                try:
                    if (tokens[0]=='float'):
                        self.params[tokens[1]] = float(tokens[2]) 
    
                    if (tokens[0]=='string'):
                        self.params[tokens[1]] = tokens[2] 
    
                    if (tokens[0]=='int'):
                        self.params[tokens[1]] = int(tokens[2])
                         
                    if (tokens[0]=='stringlist'):
                        self.params[tokens[1]] = [ token for token in tokens[2:] ]
    
                    if (tokens[0]=='stringlistbuild'):
                        self.params[tokens[1]].append(tuple(token for token in tokens[2:]))
                    
                    if (tokens[0]=='intlist'):
                        self.params[tokens[1]] = [ int(token) for token in tokens[2:] ]
    
                    if (tokens[0]=='floatlist'):
                        # adds the float list as a generator expression to the list of other generator expressions of type tokens[1]
                        # i.e. floatList CRITICALPOINT x0 y0 z0 a b c 
                        # causes a generator expression to be added to the list of generator expressions for CRITICCALPOINT in the dictionary. 
                        # each generator expression in the list of minima returns X0,y0,z0, a b and then c when it is evaluated.
                        try: 
                            self.params[tokens[1]].append( tuple(float(token) for token in tokens[2:]) )
                        except KeyError:
                            self.params[tokens[1]] = []
                            self.params[tokens[1]].append( tuple(float(token) for token in tokens[2:]) )
                except IndexError:
                    print "Ignoring Stray Line With Nothing In It"
                         
        return noErrors
    
    def dumpParams(self):
        print("Parameters were:")
        
        for key in self.__dict__:
            if ( type(self.__dict__[key]) not in [type(['a']), type({})] ):
                print(key, self.__dict__[key])
            else:
                print(key, "list")

        if  hasattr(self, 'structuralRange'):
            print("\nStructural Parameters were:\n")
            for key in self.structuralRange:
                print(key, self.structuralRange[key])
                
        return
    
    
