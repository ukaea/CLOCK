execute_process(COMMAND executable 
  TIMEOUT 1000 # it should be less than in add_test
  RESULT_VARIABLE status
)
file(REMOVE test_grey.png)
file(REMOVE test_rgb.png)