CXX=gcc

CXXFLAGS = -std=c11 -O3 -lm -Wunused-result
LIBS = -lm

OBJS = auxfunctions.o main.o initmethods.o lloyd.o elkan.o hamerly.o macqueen.o hartiganwong.o metrics.o
APP = clbin


build: $(APP)


obj/main.o: main.c
	$(CXX) $(CXXFLAGS) -c main.c -o main.o


$(OBJ)%.o: %.c
	$(CXX) $(CXXFLAGS) $(LIBS) -c $< -o $@
	
	
$(APP): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LIBS) $(OBJS) -o $(APP) 

.PHONY : clean
clean:
	rm -f $(OBJS) $(APP).exe $(APP)