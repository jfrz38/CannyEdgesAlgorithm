#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <ctime> 

#include <C_General.hpp>
#include <C_Trace.hpp>
#include <C_File.hpp>
#include <C_Arguments.hpp>
#include <C_Matrix.hpp>
#include <C_Image.hpp>

//void programa(string entrada, string salida);
void programa(char entrada[], char salida[], int tipo);
C_Matrix convolucion(C_Matrix m1, C_Matrix m2);
int direccion_cercana(C_Matrix::ElementT f);
void crear_matriz_nomax_orientacion(int i, int j);
void crear_matriz_nomax(int i, int j);
void seguir_cadena_m2(int i, int j);
void seguir_cadena_orientacion(int i, int j);

int tam,metodo,mascara;
int u_max, u_min;
int total_bien = 0;
int total_mal = 0;
C_Image imagen;
C_Matrix kernel_gauss;
C_Matrix matriz_J;
C_Matrix matriz_es;
C_Matrix matriz_eo;
C_Matrix gradiente_sobel_Jx;
C_Matrix gradiente_sobel_Jy;
C_Matrix gradiente_prewitt_Jx;
C_Matrix gradiente_prewitt_Jy;
C_Matrix gradiente_roberts_Jx;
C_Matrix gradiente_roberts_Jy;
C_Matrix mascara_filtro_x;
C_Matrix mascara_filtro_y;
C_Matrix matriz_Jx;
C_Matrix matriz_Jy;
C_Matrix matriz_direccion;
C_Matrix matriz_nomax;
C_Matrix matriz_visitados;
C_Matrix matriz_umbral;
unsigned t0, t1;
C_Image imagen_final;

int main(int argc, char **argv)
{
	printf("Inicio del programa\n");
	char txt_entrada[100];
	char txt_salida[100];
	cout << "Enter image input name:\n";
	cin.getline(txt_entrada, 100);
	cout << "Enter image output name:\n";
	cin.getline(txt_salida, 100);
	if (strcmp(txt_entrada, txt_salida) == 0)printf("Image will be overwrite\n");
	printf("Convolution kernel size (odd value and smaller than image): ");
	scanf_s("%d", &tam);
	while (true) {
		printf("Threshold max value: ");
		scanf_s("%d", &u_max);
		printf("Threshold min value: ");
		scanf_s("%d", &u_min);
		if (u_max > u_min) break;
		else printf("max value must be higher\n");
	}
	while (true) {
		cout << "Seleccionar máscara: 1 = Sobel ; 2 = Prewitt ; 3 = Roberts \n";
		scanf_s("%d", &mascara);
		if (mascara == 1 || mascara == 2 || mascara == 3) break;
		else printf("Number 1, 2 or 3\n");
	}
	while(true){
		cout << "Seleccionar rutina: 1 = vecino más cercano ; 2 = comprobar 8 vecinos:\n";
		scanf_s("%d", &metodo);
		if (metodo == 1 || metodo == 2) break;
		else printf("Number 1 or 2\n");
	}
	
	//Medir el tiempo
	t0 = clock();
	//Inicialización del programa
	programa(txt_entrada, txt_salida, metodo);

	printf("Finish successfully\n");
	
	//Finalización medición del tiempo
	t1 = clock();

	double time = (double(t1 - t0) / CLOCKS_PER_SEC);
	cout << "Execution Time: " << time << endl;
	
}

/*	Programa principal que realiza el procesamiento de la imagen
	entrada[] : nombre de la imagen a procesar
	salida[] : nombre de la imagen creada
*/
void programa(char entrada[], char salida[], int tipo) {

	//Leer imagen de entrada
	printf("Leyendo imagen\n");

	imagen.ReadBMP(entrada);
	if (imagen.Fail()) {
		cout << "Image doesn't exist.\n";
		cin.getline(entrada, 100);
		system("pause");
		exit(-1);
	}

	//Inicializar matrices
	kernel_gauss.Resize(1, tam, 1, tam);
	matriz_J.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	matriz_es.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	matriz_eo.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	gradiente_sobel_Jx.Resize(0, 2, 0, 2);
	gradiente_sobel_Jy.Resize(0, 2, 0, 2);
	gradiente_prewitt_Jx.Resize(0, 2, 0, 2);
	gradiente_prewitt_Jy.Resize(0, 2, 0, 2);
	gradiente_roberts_Jx.Resize(0, 2, 0, 2);
	gradiente_roberts_Jy.Resize(0, 2, 0, 2);
	matriz_Jx.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	matriz_Jy.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	matriz_direccion.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol());
	matriz_nomax.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	matriz_visitados.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol());
	matriz_umbral.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	

	//Calcular kernel gaussiano
	//Para ello debemos crear una matriz kernel según una variable de desviación por definir.
	//El kernel tendrá un tamaño de 5x5
	//La fórmula será: G(x,y) = (1/(2*PI*(sigma^2)))*e^-((x^2+y^2)/2+sigma^2)
	//x = fila ; y = columna ; sigma = desviación estándar
	double sigma = 1.5;
	for (int i = kernel_gauss.FirstRow(); i <= kernel_gauss.LastRow(); i++) {
		for (int j = kernel_gauss.FirstCol(); j <= kernel_gauss.LastCol(); j++) {
			kernel_gauss(i, j) = (1 / (2 * M_PI)) * exp(-((pow(i, 2) + pow(j, 2)) / (2 * pow(sigma, 2))));
		}
	}
	
	kernel_gauss.DivideEscalar(kernel_gauss.Sum());

	//Matriz_J = I * G siendo I la original y G la gaussiana creada anteriormente
	//Hacer convolución de I usando la máscara G y guardarlo en J

	printf("Convolución sobre la imagen con una máscara Gaussiana\n");
	matriz_J = convolucion(imagen, kernel_gauss);

	//Gradiente Jx y Jy

	/*
	Máscara 1 = Sobel
	Máscara 2 = Prewitt
	Máscara 3 = Roberts
	Default = Sobel
	*/
	if (mascara == 1) {
		//Máscara de Sobel
		//Gradiente X
		gradiente_sobel_Jx(0, 0) = -1;		gradiente_sobel_Jx(0, 1) = 0;		gradiente_sobel_Jx(0, 2) = 1;
		gradiente_sobel_Jx(1, 0) = -2;		gradiente_sobel_Jx(1, 1) = 0;		gradiente_sobel_Jx(1, 2) = 2;
		gradiente_sobel_Jx(2, 0) = -1;		gradiente_sobel_Jx(2, 1) = 0;		gradiente_sobel_Jx(2, 2) = 1;


		//Gradiente Y
		gradiente_sobel_Jy(0, 0) = -1;		gradiente_sobel_Jy(0, 1) = -2;		gradiente_sobel_Jy(0, 2) = -1;
		gradiente_sobel_Jy(1, 0) = 0;		gradiente_sobel_Jy(1, 1) = 0;		gradiente_sobel_Jy(1, 2) = 0;
		gradiente_sobel_Jy(2, 0) = 1;		gradiente_sobel_Jy(2, 1) = 2;		gradiente_sobel_Jy(2, 2) = 1;

		mascara_filtro_x = gradiente_sobel_Jx;
		mascara_filtro_y = gradiente_sobel_Jy;

	}
	else if (mascara == 2) {
		//Máscara de Prewitt
		//Gradiente X
		gradiente_prewitt_Jx(0, 0) = -1;		gradiente_prewitt_Jx(0, 1) = 0;		gradiente_prewitt_Jx(0, 2) = 1;
		gradiente_prewitt_Jx(1, 0) = -1;		gradiente_prewitt_Jx(1, 1) = 0;		gradiente_prewitt_Jx(1, 2) = 1;
		gradiente_prewitt_Jx(2, 0) = -1;		gradiente_prewitt_Jx(2, 1) = 0;		gradiente_prewitt_Jx(2, 2) = 1;


		//Gradiente Y
		gradiente_prewitt_Jy(0, 0) = 1;			gradiente_prewitt_Jy(0, 1) = 1;		gradiente_prewitt_Jy(0, 2) = 1;
		gradiente_prewitt_Jy(1, 0) = 0;			gradiente_prewitt_Jy(1, 1) = 0;		gradiente_prewitt_Jy(1, 2) = 0;
		gradiente_prewitt_Jy(2, 0) = -1;		gradiente_prewitt_Jy(2, 1) = -1;	gradiente_prewitt_Jy(2, 2) = -1;

		mascara_filtro_x = gradiente_prewitt_Jx;
		mascara_filtro_y = gradiente_prewitt_Jy;
	}
	else if (mascara == 3) {
		//Máscara de Roberts
		//Gradiente X
		gradiente_roberts_Jx(0, 0) = -1;	gradiente_roberts_Jx(0, 1) = 0;		gradiente_roberts_Jx(0, 2) = 0;
		gradiente_roberts_Jx(1, 0) = 0;		gradiente_roberts_Jx(1, 1) = 1;		gradiente_roberts_Jx(1, 2) = 0;
		gradiente_roberts_Jx(2, 0) = 0;		gradiente_roberts_Jx(2, 1) = 0;		gradiente_roberts_Jx(2, 2) = 0;


		//Gradiente Y
		gradiente_roberts_Jy(0, 0) = 0;		gradiente_roberts_Jy(0, 1) = 0;		gradiente_roberts_Jy(0, 2) = -1;
		gradiente_roberts_Jy(1, 0) = 0;		gradiente_roberts_Jy(1, 1) = 1;		gradiente_roberts_Jy(1, 2) = 0;
		gradiente_roberts_Jy(2, 0) = 0;		gradiente_roberts_Jy(2, 1) = 0;		gradiente_roberts_Jy(2, 2) = 0;

		mascara_filtro_x = gradiente_roberts_Jx;
		mascara_filtro_y = gradiente_roberts_Jy;
	}
	else {
		//No debería entrar aquí
		//Se aplica Sobel por defecto
		//Máscara de Sobel
		//Gradiente X
		gradiente_sobel_Jx(0, 0) = -1;		gradiente_sobel_Jx(0, 1) = 0;		gradiente_sobel_Jx(0, 2) = 1;
		gradiente_sobel_Jx(1, 0) = -2;		gradiente_sobel_Jx(1, 1) = 0;		gradiente_sobel_Jx(1, 2) = 2;
		gradiente_sobel_Jx(2, 0) = -1;		gradiente_sobel_Jx(2, 1) = 0;		gradiente_sobel_Jx(2, 2) = 1;


		//Gradiente Y
		gradiente_sobel_Jy(0, 0) = -1;		gradiente_sobel_Jy(0, 1) = -2;		gradiente_sobel_Jy(0, 2) = -1;
		gradiente_sobel_Jy(1, 0) = 0;		gradiente_sobel_Jy(1, 1) = 0;		gradiente_sobel_Jy(1, 2) = 0;
		gradiente_sobel_Jy(2, 0) = 1;		gradiente_sobel_Jy(2, 1) = 2;		gradiente_sobel_Jy(2, 2) = 1;

		mascara_filtro_x = gradiente_sobel_Jx;
		mascara_filtro_y = gradiente_sobel_Jy;

	}

	//Convolución Jx
	//Hacer convolución de J usando la máscara según el filtro elegido y guardarlo en matriz_Jx

	printf("Convolución con el gradiente X\n");
	matriz_Jx = convolucion(matriz_J, mascara_filtro_x);

	//Convolución Jy
	//Hacer convolución de J usando la máscara según el filtro elegido y guardarlo en matriz_Jy
	printf("Convolución con el gradiente Y\n");
	matriz_Jy = convolucion(matriz_J, mascara_filtro_y);

	printf("Calculando magnitud de los bordes\n");
	//Calcular la magnitud de los bordes 
	for (int i = matriz_es.FirstRow(); i <= matriz_es.LastRow(); i++) {
		for (int j = matriz_es.FirstCol(); j <= matriz_es.LastCol(); j++) {
			matriz_es(i, j) = sqrt((pow(matriz_Jx(i, j), 2) + (pow(matriz_Jy(i, j), 2))));
		}
	}

	///Método1
	if (metodo == 1) {
		
		printf("Método 1\n");
		printf("Estimando orientación de los bordes\n");
		//Estimar la orientación de la normal de los bordes
		for (int i = matriz_eo.FirstRow(); i <= matriz_eo.LastRow(); i++) {
			for (int j = matriz_eo.FirstCol(); j <= matriz_eo.LastCol(); j++) {
				matriz_eo(i, j) = atan(matriz_Jy(i, j) / matriz_Jx(i, j));
			}
		}

		//Matriz de dirección
		//Estimar dirección posible entre 0 - 45 - 90 - 135
		printf("Estimando ángulo de los bordes\n");
		for (int i = matriz_eo.FirstRow(); i <= matriz_eo.LastRow(); i++) {
			for (int j = matriz_eo.FirstCol(); j <= matriz_eo.LastCol(); j++) {
				matriz_direccion(i, j) = direccion_cercana(matriz_eo(i, j));
			}
		}

		//Matriz no_max vecino más cercano
		printf("Calculando matriz de no máximos\n");
		matriz_nomax.SetValue(0);

		for (int i = matriz_es.FirstRow(); i <= matriz_es.LastRow(); i++) {
			for (int j = matriz_es.FirstCol(); j <= matriz_es.LastCol(); j++) {
				crear_matriz_nomax_orientacion(i, j);
			}
		}

		//Histéresis según el umbral
		printf("Recorrer imagen según añadiendo valores según el umbral");
		matriz_visitados.SetValue(0);
		matriz_umbral.SetValue(0);
		for (int i = imagen.FirstRow(); i <= imagen.LastRow(); i++) {
			for (int j = imagen.FirstCol(); j <= imagen.LastCol(); j++) {
				//Si se ha visitado el punto continua la ejecución
				if (matriz_visitados(i, j) == 1)continue;
				//Los bordes se dejan igual para evitar problemas al visitar vecinos que estén fuera de rango
				if (i != imagen.FirstRow() || i != imagen.LastRow() || j != imagen.LastCol() || j != imagen.FirstCol()) {
					if (matriz_nomax(i, j) >= u_max) seguir_cadena_orientacion(i, j);
				}
			}
		}
	}
	///Método2
	else if (metodo == 2) {

		//Comprobar el peso con los vecino en la dirección correcta
		printf("Calculando matriz de no máximos\n");
		matriz_nomax.SetValue(0);

		for (int i = matriz_es.FirstRow(); i <= matriz_es.LastRow(); i++) {
			for (int j = matriz_es.FirstCol(); j <= matriz_es.LastCol(); j++) {
				crear_matriz_nomax(i, j);
			}
		}

		//Umbral histéresis

		//Histéresis
		printf("Histéresis una vez definidos el umbral\n");
		matriz_visitados.SetValue(0);
		matriz_umbral.SetValue(0);
		for (int i = imagen.FirstRow(); i <= imagen.LastRow(); i++) {
			for (int j = imagen.FirstCol(); j <= imagen.LastCol(); j++) {
				//Si se ha visitado el punto continua la ejecución
				if (matriz_visitados(i, j) == 1)continue;
				//Los bordes se dejan igual para evitar problemas al visitar vecinos que estén fuera de rango
				if (i != imagen.FirstRow() || i != imagen.LastRow() || j != imagen.LastCol() || j != imagen.FirstCol()) {
					if (matriz_nomax(i, j) >= u_max)seguir_cadena_m2(i, j);
				}
			}
		}
	}
	///Error
	else {
		printf("Error\n");
	}

	//Crear imagen final
	printf("Creando imagen final\n");
	if (matriz_umbral.Max() > 255.0) matriz_umbral.Stretch(0, 255);
	imagen_final = C_Image(matriz_umbral);
	imagen_final.palette = imagen.palette;
	imagen_final.WriteBMP(salida);

}
/*	Método para realizar convolución sobre una matriz
	m1 : matriz sobre la que se realiza la convolución
	m2 : máscara kernel
	return aux : Matriz tras la convolución
*/
C_Matrix convolucion(C_Matrix m1, C_Matrix m2) {

	int center = m2.RowN() / 2;
	C_Matrix aux(m1.FirstRow(), m1.LastRow(), m1.FirstCol(), m1.LastCol(), 255);

	C_Matrix::ElementT sumatoria_convolucion;
	for (int i = m1.FirstRow() + center; i < m1.LastRow() - center; i++) {
		for (int j = m1.FirstCol() + center; j < m1.LastCol() - center; j++) {
			sumatoria_convolucion = 0;
			for (int k = m2.FirstRow(); k <= m2.LastRow(); k++) {
				for (int l = m2.FirstCol(); l <= m2.LastCol(); l++) {
					//Se resta center por el tamaño del kernel
					sumatoria_convolucion += m1(i + k - center, j + l - center)*m2(k, l);
				}
			}
			aux(i, j) = sumatoria_convolucion;
		}
	}

	return aux;
}

/*	Método para comprobar a qué ángulo se acerca más de 0º, 42º, 90º o 135º
	f : valor a estudiar
	return : Ángulo al que más se aproxima la entrada
*/
int direccion_cercana(C_Matrix::ElementT f) {

	C_Matrix::ElementT angulo = (f / M_PI) * 180.0;
	if ((angulo < 22.5 && angulo > -22.5) || (angulo > 157.5 || angulo < -157.5)) return 0;
	if ((angulo > 22.5 && angulo < 67.5) || (angulo < -112.5 && angulo > -157.5)) return 45;
	if ((angulo > 67.5 && angulo < 112.5) || (angulo < -67.5 && angulo > -112.5)) return 90;
	if ((angulo > 112.5 && angulo < 157.5) || (angulo < -22.5 && angulo > -67.5)) return 135;
	
	return -1;	//No llega aquí
}

/*	Método para suprimir máximos sin necesidad de mirar el vecino cercano.
	i : posición i (filas) de la matriz
	j : posición j (columnas) de la matriz
*/
void crear_matriz_nomax(int i, int j) {
	
	C_Matrix::ElementT v1, v2, v3, v4;
	C_Matrix::ElementT peso, aux, aux1, aux2;

	//Evitar index out of bounds
	if (i == matriz_es.FirstRow() || i == matriz_es.LastRow() || j == matriz_es.FirstCol() || j == matriz_es.LastCol()) {
		matriz_nomax(i, j) = 0;
		return;
	}

	if (matriz_es(i, j) == 0) matriz_nomax(i, j) = 0;
	else {
		aux = matriz_es(i, j);
		//Si Jy > Jx es que tiende hacia el gradiente y
		if (abs(matriz_Jy(i, j)) > abs(matriz_Jx(i, j))) {
			if (matriz_Jy(i, j) == 0) {
				peso = 1;
			}
			else {
				peso = fabs(matriz_Jx(i, j)) / fabs(matriz_Jy(i, j));
			}
			//Variables de arriba y abajo
			v2 = matriz_es(i - 1, j);
			v4 = matriz_es(i + 1, j);

			if (matriz_Jx(i, j)*matriz_Jy(i, j)>0) {
				v1 = matriz_es(i - 1, j - 1);
				v3 = matriz_es(i + 1, j + 1);
			}
			else {
				v1 = matriz_es(i - 1, j + 1);
				v3 = matriz_es(i + 1, j - 1);
			}
		}
		else {	//Tiende hacia X
			if (matriz_Jx(i, j) == 0) {
				peso = 1;
			}
			else {
				peso = fabs(matriz_Jy(i, j)) / fabs(matriz_Jx(i, j));
			}
			v2 = matriz_es(i, j + 1);
			v4 = matriz_es(i, j - 1);

			if (matriz_Jx(i, j)*matriz_Jy(i, j) > 0) {
				v1 = matriz_es(i + 1, j + 1);
				v3 = matriz_es(i - 1, j - 1);
			}
			else {
				v1 = matriz_es(i - 1, j + 1);
				v3 = matriz_es(i + 1, j - 1);
			}
		}
		//Calculamos la interpolacion de los gradientes
		aux1 = peso*v1 + (1 - peso)*v2;
		aux2 = peso*v3 + (1 - peso)*v4;

		//Suponemos el pixel actual temp como el maximo local
		//y comprobamos si puede pertenecer al borde
		if (aux >= aux1 && aux >= aux2) {
			matriz_nomax(i, j) = aux;
		}
		else if (aux<aux1 || aux<aux2) {
			matriz_nomax(i, j) = 0;
		}
	}
}

/*	Método para crear la imágen según el umbral de mínimos de forma recursiva siguiendo los 8 vecinos
	i : posición fila de la matriz
	j : posición columna de la matriz
*/
void seguir_cadena_m2(int i, int j) {

	C_Matrix::ElementT valor;

	if(matriz_visitados(i, j) == 0)matriz_visitados(i, j) = 1;
	else return;

	matriz_umbral(i, j) = 255;

	//Vecinos del punto X
	int pos_x[8] = { 1,1,0,-1,-1,-1,0,1 };
	int pos_y[8] = { 0,1,1,1,0,-1,-1,-1 };

	//Se recorren los vecinos sin importar la dirección del punto
	for (int k = 0; k < 8; k++) {

		valor = matriz_nomax(i + pos_x[k], j + pos_y[k]);

		if (valor >= u_min) {
			seguir_cadena_m2(i + pos_x[k], j + pos_y[k]);
		}
	}
}

/*	Método para suprimir máximos dependiendo de la dirección
	i : posición i (filas) de la matriz
	j : posición j (columnas) de la matriz
*/
void crear_matriz_nomax_orientacion(int i, int j) {

	if (i == matriz_es.FirstRow() || i == matriz_es.LastRow() || j == matriz_es.FirstCol() || j == matriz_es.LastCol()) {
		matriz_nomax(i, j) = 0;
		return;
	}

	int direccion = matriz_direccion(i, j);
	switch (direccion) {
	case 0:	//Comprobar con los píxeles de la izquierda y la derecha
		if (matriz_es(i, j) < matriz_es(i, j - 1) || matriz_es(i, j) < matriz_es(i, j + 1)) {
			matriz_nomax(i, j) = 0;
		}
		else {
			matriz_nomax(i, j) = matriz_es(i, j);
		}
		break;
	case 45://Comprobar con los píxeles de la izquierda abajo y la derecha arriba
		if (matriz_es(i, j) < matriz_es(i - 1, j + 1) || matriz_es(i, j) < matriz_es(i + 1, j - 1)) {
			matriz_nomax(i, j) = 0;
		}
		else {
			matriz_nomax(i, j) = matriz_es(i, j);
		}
		break;
	case 90://Comprobar con los píxeles de arriba y abajo
		if (matriz_es(i, j) < matriz_es(i - 1, j) || matriz_es(i, j) < matriz_es(i + 1, j)) {
			matriz_nomax(i, j) = 0;
		}
		else {
			matriz_nomax(i, j) = matriz_es(i, j);
		}
		break;
	case 135://Comprobar con los píxeles de la izquierda arriba y la derecha abajo
		if (matriz_es(i, j) < matriz_es(i - 1, j - 1) || matriz_es(i, j) < matriz_es(i + 1, j + 1)) {
			matriz_nomax(i, j) = 0;
		}
		else {
			matriz_nomax(i, j) = matriz_es(i, j);
		}
		break;
	default:
		matriz_nomax(i, j) = 0;
		break;
	}
}

/*	Método para crear la imágen según el umbral de mínimos de forma recursiva según la orientación
i : posición fila de la matriz
j : posición columna de la matriz
*/
void seguir_cadena_orientacion(int i, int j) {

	if (matriz_visitados(i, j) == 1)return;
	if (i == imagen.FirstRow() || i == imagen.LastRow() || j == imagen.FirstCol() || j == imagen.LastCol()) return;

	matriz_visitados(i, j) == 1;	//Visitado
	matriz_umbral(i, j) = 255;		//Marcado como borde

	//A partir de aquí recorrer píxeles conectados en ambas direcciones perpendiculares a la normal del borde mientras sea > u_min

	//Valores de los dos vecinos
	int aux_x1, aux_y1, aux_x2, aux_y2;

	int direccion = matriz_direccion(i, j);
	switch (direccion) {
	case 0:	//Comprobar con los píxeles de la izquierda y la derecha
		aux_x1 = 0; aux_x2 = 0; aux_y1 = -1; aux_y2 = 1;
		break;
	case 45://Comprobar con los píxeles de la izquierda abajo y la derecha arriba
		aux_x1 = -1; aux_x2 = 1; aux_y1 = -1; aux_y2 = -1;
		break;
	case 90://Comprobar con los píxeles de arriba y abajo
		aux_x1 = -1; aux_x2 = 1; aux_y1 = 0; aux_y2 = 0;
		break;
	case 135://Comprobar con los píxeles de la izquierda arriba y la derecha abajo
		aux_x1 = -1; aux_x2 = 1; aux_y1 = -1; aux_y2 = 1;
		break;
	default:
		printf("Error\n");
		aux_x1 = 0; aux_x2 = 0; aux_y1 = 0; aux_y2 = 0;
		break;
	}
	//Seguir cadena por los puntos donde el valor sea mayor al umbral mínimo
	if (matriz_nomax(i + aux_x1, j + aux_y1) >= u_min) seguir_cadena_orientacion(i + aux_x1, j + aux_y1);
	if (matriz_nomax(i + aux_x2, j + aux_y2) >= u_min) seguir_cadena_orientacion(i + aux_x2, j + aux_y2);
}