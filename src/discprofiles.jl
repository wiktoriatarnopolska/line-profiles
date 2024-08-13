using Gradus, Plots
using LaTeXStrings

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
    p=plot(r,[dx_ct],fillrange=dx_cb,label=L"\dot{M} / \dot{M}_{\textrm{Edd}} = 0.3", color = :darkorchid1 ,fillcolor = :darkorchid1, xlabel=L"\textrm{Radius~} (GM/c^2)",ylabel=L"\textrm{Height~} (GM/c^2)")
    plot!(r,dx_bt,fillrange=dx_bb,label=L"\dot{M} / \dot{M}_\textrm{Edd} = 0.2", color = :hotpink2, fillcolor = :hotpink2)
    plot!(r,dx_at,fillrange=dx_ab, label=L"\dot{M} / \dot{M}_\textrm{Edd} = 0.1", color = :tan1, fillcolor = :tan1)
    x=zeros(100)
    plot!(r,x,label=L"\dot{M} / \dot{M}_\textrm{Edd} = 0.0",color=:white, aspect_ratio = 3)
    return p
end 
#a=0.0 α13=0.35
m1_a=JohannsenMetric(M=1.0,a=0.0,α13=0.35,ϵ3=0.0)
d1_at,d1_ab,d1_bt,d1_bb,d1_ct,d1_cb=disk_height(m1_a)
p1=plot_data(d1_at,d1_ab,d1_bt,d1_bb,d1_ct,d1_cb)
plot_horizon!(m1_a, color = :black, label = L"\textrm{Horizon}")

#a=0.0 α13=0.0
m2_a=JohannsenMetric(M=1.0,a=0.0,α13=0.0,ϵ3=0.0)
d2_at,d2_ab,d2_bt,d2_bb,d2_ct,d2_cb=disk_height(m2_a)
p2=plot_data(d2_at,d2_ab,d2_bt,d2_bb,d2_ct,d2_cb)
plot_horizon!(m1_a, color = :black, label = L"\textrm{Horizon}")


#a=0.0 α13=-0.35
m3_a=JohannsenMetric(M=1.0,a=0.0,α13=-0.35,ϵ3=0.0)
d3_at,d3_ab,d3_bt,d3_bb,d3_ct,d3_cb=disk_height(m3_a)

p3=plot_data(d3_at,d3_ab,d3_bt,d3_bb,d3_ct,d3_cb)
plot_horizon!(m1_a, color = :black, label = L"\textrm{Horizon}")
#a=0.8 α13=0.35

m4_a=JohannsenMetric(M=1.0,a=0.8,α13=0.35,ϵ3=0.0)
d4_at,d4_ab,d4_bt,d4_bb,d4_ct,d4_cb=disk_height(m4_a)
p4=plot_data(d4_at,d4_ab,d4_bt,d4_bb,d4_ct,d4_cb)
plot_horizon!(m1_a, color = :black, label = L"\textrm{Horizon}")
#a=0.0 α13=0.0
m5_a=JohannsenMetric(M=1.0,a=0.8,α13=0.0,ϵ3=0.0)
d5_at,d5_ab,d5_bt,d5_bb,d5_ct,d5_cb=disk_height(m5_a)
p5=plot_data(d5_at,d5_ab,d5_bt,d5_bb,d5_ct,d5_cb)
plot_horizon!(m1_a, color = :black, label = L"\textrm{Horizon}")
#a=0.0 α13=-0.35
m6_a=JohannsenMetric(M=1.0,a=0.8,α13=-0.35,ϵ3=0.0)
d6_at,d6_ab,d6_bt,d6_bb,d6_ct,d6_cb=disk_height(m6_a)
p6=plot_data(d6_at,d6_ab,d6_bt,d6_bb,d6_ct,d6_cb)
plot_horizon!(m1_a, color = :black, label = L"\textrm{Horizon}")

#a=0.998 α13=0.35

m7_a=JohannsenMetric(M=1.0,a=0.998,α13=0.35,ϵ3=0.0)
d7_at,d7_ab,d7_bt,d7_bb,d7_ct,d7_cb=disk_height(m7_a)
p7=plot_data(d7_at,d7_ab,d7_bt,d7_bb,d7_ct,d7_cb)
plot_horizon!(m1_a, color = :black, label = L"\textrm{Horizon}")
#a=0.998 α13=0.0
m8_a=JohannsenMetric(M=1.0,a=0.998,α13=0.0,ϵ3=0.0)
d8_at,d8_ab,d8_bt,d8_bb,d8_ct,d8_cb=disk_height(m8_a)
p8=plot_data(d8_at,d8_ab,d8_bt,d8_bb,d8_ct,d8_cb)
plot_horizon!(m1_a, color = :black, label = L"\textrm{Horizon}")
#a=0.998 α13=-0.35
m9_a=JohannsenMetric(M=1.0,a=0.998,α13=-0.35,ϵ3=0.0)
d9_at,d9_ab,d9_bt,d9_bb,d9_ct,d9_cb=disk_height(m9_a)
p9=plot_data(d9_at,d9_ab,d9_bt,d9_bb,d9_ct,d9_cb)
plot_horizon!(m1_a, color = :black, label = L"\textrm{Horizon}")
plot(p1,p2,p3,p4,p5,p6,p7,p8,p9,layout=@layout[p1 p2 p3; p4 p5 p6; p7 p8 p9],title=["a=0.0 α13=0.35" "a=0.0 α13=0.0" "a=0.0 α13=-0.35" "a=0.8 α13=0.35" "a=0.8 α13=0.0" "a=0.8 α13=-0.35" "a=0.998 α13=0.35" "a=0.998 α13=0.0" "a=0.998 α13=-0.35"],size=(1000,1000))