all: algorithm_plain.cpp
	g++ -std=c++11 algorithm_plain.cpp -lboost_program_options -o ap
	g++ -std=c++11 algorithm_with_NA.cpp -lboost_program_options -o ap_NA	
clean:
	$(RM) ap
	$(RM) ap_NA