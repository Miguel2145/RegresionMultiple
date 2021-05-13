/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package examples.bookTrading2;

/**
 *
 * @author migue
 */
public class MultiR {
    public static void RegresionLinealM() {
        double[] x1 = new double[] {41.9,43.4,43.9,44.5,47.3,47.5,47.9,50.2,52.8,53.2,56.7,57.0,63.5,65.3,71.1,77.0,77.8};
		double[] x2 = new double[] {29.1,29.3,29.5,29.7,29.9,30.3,30.5,30.7,30.8,30.9,31.5,31.7,31.9,32.0,32.1,32.5,32.9};
		double[] y = new double[] {251.3,251.3,248.3,267.5,273.0,276.5,270.3,274.9,285.0,290.0,297.0,302.5,304.5,309.3,321.7,330.7,349.0};
		
		int n;
		n=y.length;
		double[][] X=new double[n][3];
		double[][] XT=new double[3][n];
		
		//Obtener X
                System.out.println("Matriz X                    Matriz Y");
		for(int i=0;i<n;i++){
			for(int j=0;j<3;j++){
				if(j==0){
				X[i][j]=1;
				}else if(j==1){
			 		X[i][j]=x1[i];
				}else {
					X[i][j]=x2[i];
				} System.out.print(X[i][j]+"\t");
                                
			};System.out.print("\t "+y[i]+"");System.out.println();
		}System.out.println();
		//Obtener X'
                System.out.println("Matriz X' o transpuesta de X\n");
		for(int i=0;i<3;i++){
			for(int j=0;j<n;j++){
				if(i==0){
					XT[i][j]=1;
				}else if(i==1){
					XT[i][j]=x1[j];
				}else{
					XT[i][j]=x2[j];
				}System.out.print(XT[i][j]+"\t");
                                
			}System.out.println();
		}System.out.println();
                
		//Multiplicacion de X'X
		double[][] R=new double[3][3];
		System.out.print("\nMatriz multiplicada de X'*X \n\n");
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				for(int k=0;k<n;k++){
					R[i][j]=R[i][j] + XT[i][k]*X[k][j];
				}
				System.out.print(String.format("%.2f", R[i][j])+"    ");
			}
			System.out.print("\n");
		}
		double[][] I=invert(R);

		//Multiplicacion X'Y
		double[] R2=new double[3];
		for(int i=0;i<3;i++){
			for(int j=0;j<n;j++){
				R2[i]=R2[i] + XT[i][j]*y[j];
                                
			}System.out.print(String.format(".%.2f",R2[i] )+"    ");
		}
		
		//Imprime inversa
		System.out.print("\n\nMatrix Inversa de X'X o (X'X)^-1 \n\n");
		for (int i=0; i<3; ++i) {
                	for (int j=0; j<3; ++j){
				System.out.print(String.format("%.2f", I[i][j])+"  \t");
			}
		System.out.print("\n");
           	}

		System.out.print("\n\nMatriz multiplicada de X'Y: \n\n");
		//Imprime X'Y
		for(int i=0;i<3;i++){
			System.out.print(R2[i]+"\n");
		}

		System.out.print("\n\nEl resultado de la Regresion multiple es: \nYhat: ");
		//Multiplicacion INV(X'X)X'Y = I*R2
		double[] MLR=new double[3];
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				MLR[i]=MLR[i] + I[i][j]*R2[j];
			}
			if(i==0)System.out.print(String.format("%.2f", MLR[i])+"\t");
			if(i==1)System.out.print(String.format("%.2f", MLR[i])+"\t");
			if(i==2)System.out.print(String.format("%.2f", MLR[i])+"\t\n");
		}
	}
    public static void result(double B0,double B1, double B2) {
		System.out.println("yhat=ß0 + ß1X1 + ß1X2 \n="+String.format("%.2f", B0)+" + "+String.format("%.2f", B1)+"x1 + "+String.format("%.2f", B2)+"x2");
	}

public static double[][] invert(double a[][]){
		int n = a.length;
            	double x[][] = new double[n][n];
            	double b[][] = new double[n][n];
            	int index[] = new int[n];
            	for (int i=0; i<n; ++i)
                	b[i][i] = 1;
		
		gaussian(a, index);

		// Update the matrix b[i][j] with the ratios stored
            	for (int i=0; i<n-1; ++i)
                	for (int j=i+1; j<n; ++j)
                    		for (int k=0; k<n; ++k)
                        		b[index[j]][k]-= a[index[j]][i]*b[index[i]][k];
		
		// Perform backward substitutions
            	for (int i=0; i<n; ++i)             {
                	x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
                	for (int j=n-2; j>=0; --j){
                    		x[j][i] = b[index[j]][i];
                    		for (int k=j+1; k<n; ++k){			
                        		x[j][i] -= a[index[j]][k]*x[k][i];
                    		}
				x[j][i] /= a[index[j]][j];
                	}
            	}
     		return x;
	}
public static void gaussian(double a[][], int index[]){
		int n = index.length;
            	double c[] = new double[n];

     		// Initialize the index
            	for (int i=0; i<n; ++i)
                	index[i] = i;
	
		// Find the rescaling factors, one from each row
            	for (int i=0; i<n; ++i) {
                	double c1 = 0;
                	for (int j=0; j<n; ++j) {
                    		double c0 = Math.abs(a[i][j]);
                    		if (c0 > c1) c1 = c0;
                	}
                	c[i] = c1;
            	}	

		// Search the pivoting element from each column
            	int k = 0;
            	for (int j=0; j<n-1; ++j) {
                	double pi1 = 0;
                	for (int i=j; i<n; ++i)  {
	            	    	double pi0 = Math.abs(a[index[i]][j]);
                    		pi0 /= c[index[i]];
                    		if (pi0 > pi1) {
                        		pi1 = pi0;
                        		k = i;
                    		}
                	}
			// Interchange rows according to the pivoting order	
			int itmp = index[j];
                	index[j] = index[k];
                	index[k] = itmp;
                	for (int i=j+1; i<n; ++i) {
                    		double pj = a[index[i]][j]/a[index[j]][j];
				// Record pivoting ratios below the diagonal
                    		a[index[i]][j] = pj;	
				// Modify other elements accordingly
                    		for (int l=j+1; l<n; ++l)
                        		a[index[i]][l] -= pj*a[index[j]][l];
                	}

		}
	}//gaussiana
}
