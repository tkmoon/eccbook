// matalloc.h -- allocate general matrices
#ifndef MATALLOC_H
#define MATALLOC_H

// It should be noted that new[] does NOT set the 
// allocated space to 0!

#define CALLOCMATRIX(var,type,rows,cols) \
   { var = new type*[rows]; \
     var[0] = new type[rows*cols]; \
     for(int i=1; i < rows; i++) { \
		var[i] = var[i-1] + cols; \
	 } \
   } 

#define FREEMATRIX(var) delete[] var[0]; delete[] var;

#define CALLOCTENSOR(name,type,nmat,rows,cols) \
   {  name =  new type**[nmat]; \
      for(int i=0; i < nmat; i++) { \
       name[i] = new type*[rows]; \
       name[i][0] = new type[rows*cols]; \
       for(int j=1; j< rows; j++) { \
 	    name[i][j] = name[i][j-1]+cols; \
       } \
    } \
   } 



#define FREETENSOR(name,nmat) { \
   for(int i = 0; i < nmat; i++) {FREEMATRIX(name[i]);} \
   delete[] name;  \
   }

#define MATDUMP(name,rows,cols) \
  { cout << #name << endl; \
 for(int i = 0; i < rows; i++) { \
   for(int j=0; j < cols; j++) { \
      cout << name[i][j] << " "; \
   }; \
   cout << endl;  \
 } \
}

#define VECDUMP(name,size) \
  { cout << #name << " "; \
    for(int i = 0; i < size; i++) { \
      cout << name[i] << " "; \
    }; \
    cout << endl;  \
  } \

#define VECDUMP2(name,size,cast) \
  { cout << #name << " "; \
    for(int i = 0; i < size; i++) { \
      cout << cast(name[i]) << " "; \
    }; \
    cout << endl;  \
  } \


#endif
