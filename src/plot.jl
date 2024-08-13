using Gradus, Plots
r_start=-50
r_end=50
n=r_end-r_start
r=LinRange(r_start,r_end,n)

function Gradus.cross_section(d, ρ)
    y=Gradus.isco(m)
    if ρ < y
        return -one(typeof(ρ))
    end
    H = (3 / 2) * inv(η) * (d.Ṁ / d.Ṁedd) * (1 - sqrt(y / ρ))
    
end


function all_height(disk)
    height=zeros(r_end)
    height_b=zeros(r_end)
    for i in 1:r_end
        height[i]=cross_section(disk,i)
        height_b[i]=0-cross_section(disk,i)
    end
    a=reverse(height)
    b=reverse(height_b)
    top=vcat(a,height)
    bottom=vcat(b,height_b) 
    return top, bottom
end


function disk_height(m)
    d1=ShakuraSunyaev(m,eddington_ratio=0.1)
    d2=ShakuraSunyaev(m,eddington_ratio=0.2)
    d3=ShakuraSunyaev(m,eddington_ratio=0.3)

    d1_t,d1_b=all_height(d1)
    d2_t,d2_b=all_height(d2)
    d3_t,d3_b=all_height(d3)
    return d1_t,d1_b,d2_t,d2_b,d3_t,d3_b
end

function plot_data(dx_at,dx_ab,dx_bt,dx_bb,dx_ct,dx_cb)
    p=plot(r,[dx_ct],fillrange=dx_cb,label="Ṁ / Ṁedd = 0.3", color = :darkorchid1 ,fillcolor = :darkorchid1, xlabel="radius",ylabel="height")
    plot!(r,dx_bt,fillrange=dx_bb,label="Ṁ / Ṁedd = 0.2", color = :hotpink2, fillcolor = :hotpink2)
    plot!(r,dx_at,fillrange=dx_ab, label="Ṁ / Ṁedd = 0.1", color = :tan1, fillcolor = :tan1)
    return p
end  
#a=0.0 α13=0.35
m1_a=JohannsenMetric(M=1.0,a=0.0,α13=0.0,ϵ3=0.0)
d1_at,d1_ab,d1_bt,d1_bb,d1_ct,d1_cb=disk_height(m1_a)
plot_data(d1_at,d1_ab,d1_bt,d1_bb,d1_ct,d1_cb)
plot_horizon!(m1_a, color = :black, label = "horizon")
