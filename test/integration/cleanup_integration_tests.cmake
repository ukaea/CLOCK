execute_process(COMMAND executable 
  TIMEOUT 1000 # it should be less than in add_test
  RESULT_VARIABLE status
)
file(REMOVE make_me_big.png)
file(REMOVE make_me_little.png)
file(REMOVE made_you_little.png)
file(REMOVE made_you_big.png)
file(REMOVE test_saltandpepperfilter.sandp1.png)
file(REMOVE test_saltandpepperfilter.sandp1_out.png)
file(REMOVE test_saltandpepperfilter.sandp3.png)
file(REMOVE test_saltandpepperfilter.sandp3_out.png)
file(REMOVE test_rgb.png)
file(REMOVE test_grey.png)