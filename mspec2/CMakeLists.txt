add_library(fabm_models_hereon OBJECT
            hereon_model_library.F90
            mspec.F90
            phy.F90
            zoo.F90
           )

add_dependencies(fabm_models_hereon fabm_base)