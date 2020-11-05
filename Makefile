CXX = g++
CXXFLAGS = 

DEF_SRC =  main.cpp matrix.cpp
DEF_T_SRC =  main.cpp matrix.cpp tests.cpp unit_tests.cpp
DEF_G_SRC = main.cpp matrix.cpp generator.cpp

work: $(DEF_SRC)
	$(CXX) $(CXXFLAGS) $? -DWORK -o prog
test: $(DEF_T_SRC)
	$(CXX) $(CXXFLAGS) $? -DTEST -o prog
generator: $(DEF_G_SRC)
	$(CXX) $(CXXFLAGS) $? -DGENERATOR -o prog
