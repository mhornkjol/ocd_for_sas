patient="PatID-004_cut"
method="ocdr"
img="stokes_vel"

if [ "${img}" = "data" ]
then
    # DATA
    declare -a dts=("0h")
    declare -a betas=(1e-2 1e-6) # betas=(1e-2 1e-4 1e-6)
    declare -a resolution=(1 4) # resolution=(1 2 4)
    declare -a dispersion=(1)
else
    # Manufactured
    declare -a dts=("0_165")
    declare -a betas=(1e-4 1e-6) #betas=(1e-4 1e-5 1e-6)
    declare -a resolution=(1 4)
    declare -a dispersion=(1)
fi

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
                    mkdir -p ${patient}/saga_pvd/${img}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}
                    scp marhorn@saga:/cluster/projects/nn9279k/Martin/sastransport_systematic/${img}/${patient}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}_1*/results*/pvd/phi.pvd ${patient}/saga_pvd/${img}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}/phi.pvd
                    scp marhorn@saga:/cluster/projects/nn9279k/Martin/sastransport_systematic/${img}/${patient}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}_1*/results*/pvd/phi000000.vtu ${patient}/saga_pvd/${img}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}/phi000000.vtu
                    scp marhorn@saga:/cluster/projects/nn9279k/Martin/sastransport_systematic/${img}/${patient}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}_1*/results*/pvd/c.pvd ${patient}/saga_pvd/${img}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}/c.pvd
                    scp marhorn@saga:/cluster/projects/nn9279k/Martin/sastransport_systematic/${img}/${patient}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}_1*/results*/pvd/c000000.vtu ${patient}/saga_pvd/${img}/pat_${patient}_n_${res}_beta_${beta}_key_${dt}_disp_${disp}_method_${method}/c000000.vtu
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
                    mkdir -p ${patient}/saga_pvd/${img}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}
                    scp marhorn@saga:/cluster/projects/nn9279k/Martin/sastransport_systematic/${img}/${patient}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}_1*/results*/pvd/phi.pvd ${patient}/saga_pvd/${img}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}/phi.pvd
                    scp marhorn@saga:/cluster/projects/nn9279k/Martin/sastransport_systematic/${img}/${patient}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}_1*/results*/pvd/phi000000.vtu ${patient}/saga_pvd/${img}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}/phi000000.vtu
                    scp marhorn@saga:/cluster/projects/nn9279k/Martin/sastransport_systematic/${img}/${patient}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}_1*/results*/pvd/c.pvd ${patient}/saga_pvd/${img}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}/c.pvd
                    scp marhorn@saga:/cluster/projects/nn9279k/Martin/sastransport_systematic/${img}/${patient}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}_1*/results*/pvd/c000000.vtu ${patient}/saga_pvd/${img}/pat_${patient}_n_${res}_beta_${beta}_dt_${dt}_disp_${disp}_method_${method}/c000000.vtu

                done
            done
        done
    done
fi