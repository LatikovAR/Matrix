CXX = g++
CXXFLAGS = 

DEF_SRC =  main.cpp matrix.cpp
DEF_T_SRC =  main.cpp matrix.cpp tests.cpp unit_tests.cpp

work: $(DEF_SRC)
	$(CXX) $(CXXFLAGS) $? -DWORK -o prog
test: $(DEF_T_SRC)
	$(CXX) $(CXXFLAGS) $? -DTEST -o prog
