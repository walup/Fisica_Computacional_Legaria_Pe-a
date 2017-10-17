__precompile__()
module herramientas

using SymPy
export newtonsMethod,newtonsMethodInterval,searchRoots,newtonsMethodSpitRoots,riemannIntegral,trapezoidIntegral,simpsonsIntegral,eulersMethod,implicitEuler,rungeKutta4,last
#=Metodo de Newton con argumentos: f, su derivada y la condicion inicial=#
function newtonsMethod(f::Function,g::Function,seed)
    MAX_N = 200;
    numIter = 0
    x = seed 
   while(abs(f(x))>0.00001)
        if(numIter<MAX_N)
        x = x-f(x)/g(x)        
        numIter+=1
        else
           break 
        end   
end 
    return x
end

#=Metodo de Newton que toma como argumento solo la funcion, de la cual se quieren obtener los ceros.=#
function newtonsMethod(f::Function,seed)
    A,x,a,n,m=symbols("A,x,a,n,m")
    derPy =simplify(diff(f),[x])
    der = lambdify(derPy,[x])
    MAX_N = 200
    valor = seed
    iteraciones = 0
    
    while(abs(f(valor))>0.00001)
        if(iteraciones<MAX_N)
            valor=valor-(f(valor)/der(valor))
            iteraciones+=1
            
        else
            println("demasiadas iteraciones, es posible que no converja")
           break 
        end
    end    
    return valor
end

#=Funcion, que aplica el metodo de Newton a todo un intervalo. Esta recibe como argumentos la funcion y su derivada. =#
function newtonsMethodInterval(f::Function,g::Function,interval)
   
    values = []
    
    for i in 1:length(interval)
       push!(values,newtonsMethod(f,g,interval[i]))
    end
    return values
end
#=Misma funcion pero solo recibe como argumento la funcion=#
function newtonsMethodInterval(f::Function,interval)

values = []
    
    for i in 1:length(interval)
       push!(values,newtonsMethod(f,interval[i]))
    end
    return values
end

function searchRoots(f::Function, epsilon, values)
    roots = []
    isEqual = false
    push!(roots,values[1])
    
    for v in values
        #=Solo vamos a considerar como candidatos a aquellas raices que evaluadas en la funcion den un numero cercano a cero=#
        
        if(abs(f(v))<0.00001)
            for v2 in roots 
                if(v2<v+epsilon && v2>v-epsilon)
                    isEqual = true                    
                end 
            end 
            
            if(isEqual)
                isEqual = false
            else
                push!(roots,v)
            end 
        end 
    end 
    return roots 
end 

function newtonsMethodSpitRoots(f::Function,g::Function,interval,epsilon)
    
    #=Para asegurarnos de que realmente nos estan dando un intervalo checamos las dimensiones de interval =#
    
    if(ndims(interval) ==1)
    values = []
    
    for i in 1:length(interval)
       push!(values,newtonsMethod(f,g,interval[i]))
    end
    
    #=Ahora tenemos que buscar las raices distintas con nuestro arreglo creado, asi que usaremos la funcion auxiliar searchRoots=#
        
        roots =  searchRoots(f,epsilon,values)
        
        return roots
        
    else
        println("Error: Se requiere dar un arreglo para esta funcion")
    end

end

function newtonsMethodSpitRoots(f::Function,interval,epsilon)
    
    #=Para asegurarnos de que realmente nos estan dando un intervalo checamos las dimensiones de interval =#
    
    if(ndims(interval) ==1)
    values = []
    
    for i in 1:length(interval)
       push!(values,newtonsMethod(f,interval[i]))
    end
    
    #=Ahora tenemos que buscar las raices distintas con nuestro arreglo creado, asi que usaremos la funcion auxiliar searchRoots=#
        
        roots =  searchRoots(f,epsilon,values)
        
        return roots
        
    else
        println("Error: Se requiere dar un arreglo para esta funcion")
    end
end  

#=Aqui comienzan los métodos de integración.=#

function riemannIntegral(f::Function,a,b)
    
    #=Esta es la implementacion del metodo de Riemann. subdividiremos el intervalo [a,b] en 200 intervalos, para obtener un valor mas exacto de la integral. Comenzamos guardando en una variable el tamaño delta de cada subintervalo. definimos las variables necesarias y comenzamos a iterar. cada ciclo del for corresponde a trabajar con uno de los 200 subintervalos. para cada uno de ellos se deben guardar los extremos del intervalo y el valor a la mitad en variables, y luego sumar a la integral la funcion evaluada en el punto medio por la longitud del intervalo delta. Una vez que se ha completado el for se regresa el valor de la integral=#
    n = 200
    delta = (b-a)/n
    xNext=0
    xBefore=0
    xMiddle=0
    
    integral = 0
    
    for i in 1:n
        xNext = a +(i)*delta
        xBefore = a +(i-1)*delta
        xMiddle = (xBefore +xNext)/2
        
        integral += f(xMiddle)*delta
    end 
    
    return integral 
end 

function trapezoidIntegral(f::Function,a,b)
    n = 200
    delta = (b-a)/n
    xNext=0
    xBefore=0
    
    integral = 0
    
    for i in 1:n
        
        xNext = a+i*delta
        xBefore = a+(i-1)*delta
        
        integral += ((f(xNext) +f(xBefore))/2)*delta
    end 
    
    return integral 
end 

function simpsonsIntegral(f::Function,a,b)
    n = 200
    delta = (b-a)/n
    xBefore=0
    xMiddle=0
    xAfter=0
    
    integral = 0
    
    for i in 1:n
        xBefore = a+(i-1)*delta
        xAfter = a+i*delta
        xMiddle = (xBefore +xAfter)/2
        
        integral += (f(xBefore) +f(xAfter) +4*f(xMiddle))*(delta/6)
    end
    return integral
end 

function eulersMethod(derivada::Function,x0,xf,y0,h)
    
    listaX = []
    listaY = []
    
    #Agreagamos los primeros elementos a las listas 
    
    push!(listaX,x0)
    push!(listaY,y0)
    
    #=Mientras no se llegue al ultimo elemento se seguiran agregando elementos a las listas. Aqui tambien se podria haber calculado el número de pasos para llegar a xf como n = (xf-x0)/h y luego implementado un for de 1 a n, pero en realidad, es lo mismo, y además si hacemos lo anterior hay que redondear a k al entero menor. i.e. debemos ocupar mas espacio para almacenar a n , y hacer mas operaciones de tiempo constante. =#
    
    while(last(listaX)<xf)
        
        push!(listaY,last(listaY) + derivada(last(listaY),last(listaX))*h)
        push!(listaX,last(listaX)+h)
    end
    
    return [listaX,listaY]
    
end

function implicitEuler(f::Function,t0::Float64,tf::Float64,x0::Float64,delta::Float64)
    #=Definimos la funcion que meteremos al método de Newton para aproximar el siguiente punto. a, se sustituira por el punto previo. =#

    a = Sym("a")
    t=Sym("t")
    g(x,t,a) = x-(f(x,t)*delta)-a
    
    
    xArray=[]
    tArray=[]
    newX = 0
    
    #=Se agregan las condiciones iniciales a los arreglos. =#
    
    push!(tArray,t0)
    push!(xArray,x0)
    
    while(last(tArray)<tf)
        
        #Se aproxima el punto siguiente ocn el método de Newton y se guarda en la variable newX
        
        newX = newtonsMethod(x->subs(g(x,t,a),t=>last(tArray),a=>last(xArray)),last(xArray))
        
        #=Una vez que se tiene el X nuevo agregamos el nuevo punto calculado con euler al arreglo x=#
        
        push!(xArray,last(xArray)+f(newX,last(tArray)+delta)*delta)
        push!(tArray,last(tArray)+delta)
    end
    
    return [tArray,xArray]
end

function rungeKutta4(f::Function,t0,tf,x0,delta)
    tArray = []
    xArray = []
    
    #=Como siempre agregamos los primeros puntos a los arreglos=#

    push!(tArray,t0)
    push!(xArray,x0)
    
    #=Definimos las funciones de Runge Kutta=#
    
    f1(x,t) = f(x,t)
    f2(x,t) = f(x+(delta/2)*f1(x,t),t+(delta/2))
    f3(x,t) = f(x+(delta/2)*f2(x,t),t+(delta/2))
    f4(x,t) = f(x+(delta)*f3(x,t),t+delta)
    
    rungeFunction(x,t) = x+(delta/6)*(f1(x,t)+2*f2(x,t) +2*f3(x,t)+f4(x,t))
    
    while(last(tArray)<tf)
        push!(xArray,rungeFunction(last(xArray),last(tArray)))
        push!(tArray,last(tArray)+delta)
        
    end
    
    return [tArray,xArray]
    
end

function last(A)
    
    return A[length(A)]
end

end
