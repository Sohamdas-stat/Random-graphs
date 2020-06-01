#Section 2.1
# 1
# this gives histogram of probability of d-neighbourhood being a tree for different d and n fixing c
library(igraph)

check.tree <- function(n, c, d){
	p = c / n
	vector = NULL
	for(i in 1:100){
		con = NULL
		for(j in 1:20){
			g = sample_gnp(n,p)
			v = sample(1:n,1)
			sub = induced_subgraph(g, which(distances(g, v, as.vector(V(g))) <= d))
			con = c(con,as.numeric(length(V(sub)) == length(E(sub)) + 1))
		}
		vector = c(vector,mean(con))
	}
	return(vector)
}

main <- function(d = 4){
	nseq = c(1000,3000,6000,10000,30000,60000)
	cseq = c(0.5,1,1.5,2)
	for (c in cseq){
		quartz(width = 10, height = 7)
		layout(matrix(1:6, nrow = 2, byrow = T))
		for (n in nseq){
			vector = check.tree(n, c, d)
			print(range(vector))
			x=vector[1:20]
			t = (mean(x)-1)/sqrt(var(x)/length(x))
			p = pt(t,length(x)-1,lower.tail=T)
			print(p)
			l = sort(vector)[3]
			u = sort(vector)[97]
			hist(vector, breaks=seq(0,1,by=.05), prob=TRUE, xlab="p_tree", col=rgb(0,0,1,1/4), main=paste("For n = ",n))
			legend("topleft", c(paste("CI = (",l,",",u,")"),paste("s.e. = ",se(vector))))
		}
	}
}

main()

#Section 2.2
# 2
# this gives the expected no. of cycles in d-neighbourhood
library(igraph)

myplot <- function(matrix){
	r = length(matrix[,1])
	c = length(matrix[1,])
	x = row.names(matrix)
	plot(NULL, xlim = c(0,10000), ylim = c(min(matrix),max(matrix)), xlab = "n", ylab = "count")
	for (i in 1:c){
		lines(x, matrix[,i], col = i, lty = 1, lwd = 2)
	}
}

find_cycle <- function(n, c, d){
	p = c / n
	g = sample_gnp(n,p)
	v = sample(1:n,1)
	sub = induced_subgraph(g, which(distances(g, v, as.vector(V(g))) <= d))
	Cycles = NULL
	for(v1 in V(sub)){
		for(v2 in neighbors(sub, v1, mode="out")){
			Cycles = c(Cycles, lapply(all_simple_paths(sub, v2,v1, mode="out"), function(p) c(v1,p)))
		}
	}
	u = length(Cycles[which(sapply(Cycles, length) == 4)])/6
	v = length(Cycles[which(sapply(Cycles, length) == 5)])/8
	w = length(Cycles[which(sapply(Cycles, length) == 6)])/10
	return(sum(u,v,w))
}

main <- function(d=4){
	nseq = c(100,300,600,1000,3000,6000)
	cseq = c(0.5,1,1.5,2)
	tmatrix = data.frame()
	for (n in nseq){
		tvector = NULL
		for (c in cseq){
			tval = 0
			for(rept in 1:500){
				v = find_cycle(n, c, d)
				tval = tval + v
			}
			tvector = c(tvector, tval/500)
		}
		tmatrix = rbind(tmatrix, tvector)
	}
	rownames(tmatrix) = nseq
	colnames(tmatrix) = cseq
	dev.new()
	myplot(tmatrix)
}

main()

# Section 2.3
# 3
# this gives histogram of probability of d-neighbourhood being a tree for different c and n fixing d
library(igraph)

check.tree <- function(n, c, d){
	p = c / n
	vector = NULL
	for(i in 1:100){
		con = NULL
		for(j in 1:20){
			g = sample_gnp(n,p)
			v = sample(1:n,1)
			sub = induced_subgraph(g, which(distances(g, v, as.vector(V(g))) <= d))
			con = c(con,as.numeric(length(V(sub)) == length(E(sub)) + 1))
		}
		vector = c(vector,mean(con))
	}
	return(vector)
}

main <- function(){
	nseq = c(2000,5000,10000,50000,100000,200000)
	dseq = seq(3,6,by=1)
	c = 3
	for (d in dseq){
		quartz(width = 10, height = 7)
		layout(matrix(1:6, nrow = 2, byrow = T))
		for (n in nseq){
			vector = check.tree(n, c, d)
			print(range(vector))
			print(mean(vector))
			l = sort(vector)[3]
			u = sort(vector)[97]
			hist(vector, breaks=seq(-0.05,1.05,by=0.1), prob=TRUE, xlab="p_tree", col=rgb(0,0,1,1/4), main=paste("For n = ",n))
			legend("topleft", c(paste("CI = (",l,",",u,")"),paste("mean = ",mean(vector))))
			legend("topright", c(paste("CI = (",l,",",u,")"),paste("mean = ",mean(vector))))
		}
	}
}

main()

#Section 2.4
# 4
# this gives histogram of degrees of randomly chosen vertices for different c and n fixing d
library(igraph)

deg <-function(n, c, d){
	vector = NULL
	p = c / n
	for (i in 1:1000){
		g = sample_gnp(n,p)
		v = sample(1:n,1)
		dg = degree(g, v)
		vector = c(vector, dg)
	}
	return(vector)
}

main <- function(d = 4){
	nseq = c(100,500,1000,10000)
	cseq = c(0.5,1,1.5,2,3,4)
	for (c in cseq){
		quartz(width = 10, height = 7)
		layout(matrix(1:6, nrow = 2, byrow = T))
		for (n in nseq){
			vector = deg(n,c,d)
			poi = round(dpois(c(1:10),c)*1000)
			poi = c(1000 - sum(poi), poi)
			poi = rep(c(0:10),poi)
			chi = 0
			for(i in 1:10){
				chi = chi + (length(which(vector == i))/1000 - dpois(i,c))^2/dpois(i,c)
			}
			p1 = round(pchisq(chi,9, lower.tail=F),2)
			x = NULL
			y = NULL
			for(i in 0:10){
				x = c(x, length(which(vector <= i))/1000)
				y = c(y, length(which(vector <= i))/1000)
			}
			p2 = round(ks.test(x,y)$p.value,2)
			p3 = round(t.test(vector,poi, alternative = "two.sided", var.equal = FALSE)$p.value,2)
			hist(vector, breaks=seq(-.5,15.5,by=1), ylim=c(0,0.1+max(dpois(1:10,c))), prob=T, xlab="degree", ylab="probability", col=rgb(0,0,1,1/4), main=paste("For n = ",n))	## range of x (-.5,15.5) for c = 4 and (-.5,12.5) for c = 3
			hist(poi, breaks=seq(-.5,15.5,by=1), ylim=c(0,max(0.1+dpois(1:10,c))), prob=T, xlab="degree", ylab="probability", col=rgb(1,0,0,1/4), add=T)
			legend("topright", c("p-value",paste("chisq.test = ",p1),paste("ks.test = ",p2),paste("t.test = ",p3)))
		}
	}
}

main()

#Section 2.5
# 5
# this gives histogram of degrees of randomly chosen vertices at t-depth in d-neighbourhood for different c and n fixing d
library(igraph)

nbdeg <- function(n, c, d, t){
	p = c / n
	vector = NULL
	for (i in 1:1000){
		g = sample_gnp(n,p)
		v = sample(1:n,1)
		ngb = which(distances(g,v,V(g)) == t)
		if(length(ngb) == 1){
			vector = c(vector, degree(g,ngb))
		}
		if(length(ngb) > 1){
			dg = degree(g, sample(ngb,1))
			vector = c(vector, dg)
		}
	}
	return(vector)
}

main <- function(d = 5){
	nseq = c(100,500,1000,10000)
	cseq = c(0.5,1,1.5,2,3,4)
	t = 1		#fix t here
	for (c in cseq){
		quartz(width = 10, height = 7)
		layout(matrix(1:4, nrow = 2, byrow = T))
		for (n in nseq){
			vector = nbdeg(n, c, d, t)
			poi = round(dpois(c(1:10),c)*1000)
			poi = c(1000 - sum(poi), poi)
			poi = rep(c(0:10),poi)
			poi = poi + 1
			chi = 0
			for(i in 1:10){
				chi = chi + (length(which(vector == i))/1000 - dpois(i,c))^2/dpois(i,c)
			}
			p1 = round(pchisq(chi,9, lower.tail=F),2)
			x = NULL
			y = NULL
			for(i in 0:10){
				x = c(x, length(which(vector <= i))/1000)
				y = c(y, length(which(vector <= i))/1000)
			}
			p2 = round(ks.test(x,y)$p.value,2)
			p3 = round(t.test(vector,poi, alternative = "two.sided", var.equal = FALSE)$p.value,2)
			hist(vector, breaks=seq(-.5,15.5,by=1), ylim=c(0,0.1+max(dpois(1:10,c))), prob=T, xlab="degree", col=rgb(0,0,1,1/4), main=paste("For n = ",n))	## range of x (-.5,15.5) for c = 4 and (-.5,12.5) for c = 3
			hist(poi, breaks=seq(-.5,15.5,by=1), ylim=c(0,max(0.1+dpois(1:10,c))), prob=T, xlab="degree", ylab="probability", col=rgb(1,0,0,1/4), add=T)
			legend("topright", c("p-value",paste("chisq.test = ",p1),paste("ks.test = ",p2),paste("t.test = ",p3)))
		}
	}
}

main()

#Section 2.6
# 6
# this gives different tests for independence of two randomly selected vertices in d-neighbourhood
library(igraph)

ngbhdeg <-function(n, c, d){
	vec1 = NULL
	vec2 = NULL
	p = c / n
	i = 0
	while(i < 100){
		g = sample_gnp(n,p)
		v = sample(1:n,1)
		dis = distances(g,v,V(g))
		if(length(which(dis <= d)) >= 5){
			vertex = sample(intersect(which(dis>0),which(dis<d)), 2, replace = FALSE)
			vec1 = c(vec1, degree(g,vertex[1]))
			vec2 = c(vec2, degree(g,vertex[2]))
			i = i + 1
		}
	}
	return(list(vec1, vec2))
}

test.kendall <- function(x,y){
	n = length(x)
	con = 0
	for(i in 1:n){
		for(j in 1:n){
			if(isTRUE((x[i]-x[j])*(y[i]-y[j]) > 0)){
				con = +1
			}
		}
	}
	r = con/(2*choose(n,2))
	return(r)
}

test.spearman <- function(x,y){
	n = length(x)
	d = rank(x)-rank(y)
	s = sum(d*d)
	r = 1-6*s/(n*(n^2-1))
	return(r)
}

test.chisq <- function(x,y){
	k = 5
	table = matrix(0,k,k)
	n = length(x)
	for(i in 1:n){
		if(x[i] < k && y[i] < k && x*y > 0){
			table[x[i],y[i]] = table[x[i],y[i]] + 1
		}
	}
	for(i in 1:(k-1)){
		table[k,i] = sum(table[,i])
		table[i,k] = sum(table[i,])
	}
	table[k,k] = sum(table[k,])
	n = table[k,k]
	chi = 0
	for(i in 1:(k-1)){
		for(j in 1:(k-1)){
			e = table[i,k]*table[k,j]/n
			chi = chi + (table[i,j] - e)^2/e
		}
	}
	p = round(pchisq(chi,(k-2)^2,lower.tail=FALSE),2)
	return(p)
}

main <- function(d = 4){
	nseq = c(seq(1000,9000,by=4000),seq(30000,90000,by=20000))
	cseq = c(1,1.5,2,2.5)
	chi = NULL
	quartz(width = 10, height = 7)
	layout(matrix(1:4, nrow = 2, byrow = T))
	for (c in cseq){
		r.k = NULL
		r.s = NULL
		p.chi = NULL
		for (n in nseq){
			out = ngbhdeg(n, c, d)
			vec1 = unlist(out[1])
			vec2 = unlist(out[2])
			r.k = c(r.k,test.kendall(vec1,vec2))
			r.s = c(r.s,test.spearman(vec1,vec2))
			p.chi = c(p.chi,test.chisq(vec1,vec2))
		}
		chi = rbind(chi,p.chi)
		print(chi)
		plot(nseq, abs(r.k), ylim = c(0,1), xlab="n", ylab="correlation", col=2, pch=20, main=paste("For c = ",c))
		points(nseq, abs(r.s), col=3, pch = 20)
		legend("topleft",legend=c("Kendall corr","Pearson rank corr"),col=c(2,3),pch=20)
		legend("topright",legend=c("p-value of chisq test"),col=c(4),pch=20)
	}
	print(chi)
	quartz(width = 10, height = 7)
	layout(matrix(1:4, nrow = 2, byrow = T))
	for(i in 1:length(chi[,1])){
		plot(nseq, chi[i,], ylim = c(0,1), xlab="n", ylab="p-value", col=4, pch=20, main=paste("For c = ",cseq[i]))
	}
}

main()

# Section 2.6
#7
# this gives the histogram of sum of degrees of two randomly selected vertices in d-neighbourhood
library(igraph)

ngbhdeg <-function(n, c, d){
	vec1 = NULL
	vec2 = NULL
	p = c / n
	i = 0
	while(i < 1000){
		g = sample_gnp(n,p)
		v = sample(1:n,1)
		dis = distances(g,v,V(g))
		if(length(which(dis <= d)) >= 5){
			vertex = sample(intersect(which(dis>0),which(dis<d)), 2, replace = FALSE)
			vec1 = c(vec1, degree(g,vertex[1]))
			vec2 = c(vec2, degree(g,vertex[2]))
			i = i + 1
		}
	}
	return(list(vec1, vec2))
}

main <- function(d = 4){
	nseq = c(500,1000,5000,10000)
	cseq = c(1,1.5,2,3)
	for (c in cseq){
		quartz(width = 10, height = 7)
		layout(matrix(1:4, nrow = 2, byrow = T))
		r.k = NULL
		r.s = NULL
		p.chi = NULL
		for (n in nseq){
			out = ngbhdeg(n, c, d)
			vec1 = unlist(out[1])
			vec2 = unlist(out[2])
			poi = round(dpois(c(1:13),2*c)*1000)
			poi = c(1000 - sum(poi), poi)
			poi = rep(c(0:13),poi)
			poi = poi + 2
			vector = vec1 + vec2
			chi = 0
			for(i in 2:15){
				chi = chi + (length(which(vector == i))/1000 - dpois(i-2,2*c))^2/dpois(i-2,2*c)
			}
			p1 = round(pchisq(chi,13, lower.tail=F),2)
			x = NULL
			y = NULL
			for(i in 0:13){
				x = c(x, length(which(vector <= i))/1000)
				y = c(y, length(which(poi <= i))/1000)
			}
			p2 = round(ks.test(x,y)$p.value,2)
			p3 = round(t.test(vector,poi, alternative = "two.sided", var.equal = FALSE)$p.value,2)
			hist(vector, breaks=seq(-.5,20.5,by=1), ylim=c(0,0.1+max(dpois(1:10,c))), prob=T, xlab="degree", col=rgb(0,0,1,1/4), main=paste("For n = ",n))	## range of x (-.5,15.5) for c = 4 and (-.5,12.5) for c = 3
			hist(poi, breaks=seq(-.5,20.5,by=1), ylim=c(0,max(0.1+dpois(1:10,c))), prob=T, xlab="degree", col=rgb(1,0,0,1/4), add=T)
			legend("topright", c("p-value",paste("chisq.test = ",p1),paste("ks.test = ",p2),paste("t.test = ",p3)))
		}
	}
}

main()
