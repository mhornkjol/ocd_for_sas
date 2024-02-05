patient="PatID-004_cut"
method="ocdr_z_vel"
img="z_vel_noise_10"

if [ "${img}" = "data" ]
then
    # DATA
    declare -a dts=("0h")
    declare -a betas=(1e-2 1e-4 1e-6)
    declare -a resolution=(1 2 4) # resolution=(1 2 4)
    declare -a dispersion=(1)
else
    # Manufactured
    declare -a dts=("0_42" "0_83" "0_125" "0_165")
    declare -a betas=(1e-4 1e-5 1e-6)
    declare -a resolution=(1 2 4) # resolution=(1 2 4)
    declare -a dispersion=(1)
fi
# declare -a dispersion=(1)

if [ "${img}" = "data" ]
then
    for dt in "${dts[@]}"
    do
        for beta in "${betas[@]}"
        do
            for disp in "${dispersion[@]}"
            do
                for res in "${resolution[@]}"
                do
                    mkdir -p ${patient}/saga_files/${img}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}
                    scp marhorn@saga:/cluster/projects/nn9279k/Martin/sastransport_systematic/${img}/${patient}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}_1*/results*/optimization_values.csv ${patient}/saga_files/${img}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}/optimization_values.csv
                done
            done
        done
    done
else
    for dt in "${dts[@]}"
    do
        for beta in "${betas[@]}"
        do
            for disp in "${dispersion[@]}"
            do
                for res in "${resolution[@]}"
                do
                    mkdir -p ${patient}/saga_files/${img}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}
                    scp marhorn@saga:/cluster/projects/nn9279k/Martin/sastransport_systematic/${img}/${patient}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}_1*/results*/optimization_values.csv ${patient}/saga_files/${img}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}/optimization_values.csv
                done
            done
        done
    done
fi