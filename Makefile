CFLAGS=-g -Wall
OBJS=main.o protein-markov-model.o protein-scanner.o dna-scanner.o
sequence-mine: $(OBJS)
	$(CXX) -o sequence-mine $(OBJS) -lpthread -lre2
%.o: %.cc
	$(CXX) $(CFLAGS) -c $<