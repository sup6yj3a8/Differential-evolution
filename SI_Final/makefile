TARGET  = main
SRC_DIR = ./
INC_DIR = ./
CC = g++

CFLAGS += -I $(INC_DIR)
SRC     = $(wildcard $(SRC_DIR)*.cpp)
OBJ     = $(SRC:%.cpp=%.o) #$(patsubst %.cpp, %.o, $(SRC))
DEP     = $(SRC:%.cpp=%.d) #$(patsubst %.o, %.d, $(OBJ))

all : $(TARGET)

$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@

# Include all .d files
-include $(DEP)

# Build target for every single object file.
# The potential dependency on header files is covered by calling `-include $(DEP)`.
%.o : %.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@

.PHONY: clean exe
# clean
clean :
	rm -f $(TARGET) $(OBJ) $(DEP)

# exe
exe :
	@./$(TARGET)
