#include <CUnit/CUnit.h>
#include "CUnit/Basic.h"
#include <iostream>
#include <vector>
#include "data.hpp"
#include "LSH.hpp"

void testDistance(void){

    vector<double> x = {10.1, 10.1, 10.1};
    vector<double> y = {10.1, 10.1, 10.1};

    int d = 3;

    CU_ASSERT(distance(x, y, d) == 0.0);

}

void testDiscrete_Frechet_Distance(void){

    vector<double> x = {5.0, 10.0, 4.0};
    vector<double> y = {3.0, 2.0, 1.0};

    CU_ASSERT(Discrete_Frechet_Distance(x, y) == 7.0);

}

void testModulo(void){
    
    LSH object = LSH(); 
    int rh = -200;
    int M = 18;

    CU_ASSERT(object.modulo(rh, M) == 16);

}

int main(){

    CU_pSuite pSuite = NULL;
   
    /* initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();
    
    /* add a suite to the registry */
    pSuite = CU_add_suite("Suite_1", NULL, NULL);
    
    if (NULL == pSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* add the tests to the suite */
    /* NOTE - ORDER IS IMPORTANT - MUST TEST fread() AFTER fprintf() */
    if ((NULL == CU_add_test(pSuite, "test of distance()", testDistance)) || (NULL == CU_add_test(pSuite, "test of modulo()", testModulo)) || (NULL == CU_add_test(pSuite, "test of Discrete_Frechet_Distance()", testDiscrete_Frechet_Distance))){
        CU_cleanup_registry();
        return CU_get_error();
    }
    
    /* Run all tests using the CUnit Basic interface */
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();

}