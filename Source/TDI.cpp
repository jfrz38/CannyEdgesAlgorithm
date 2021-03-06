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

void programa(char entrada[], char salida[]);
C_Matrix convolucion(C_Matrix m1, C_Matrix m2);
int direccion_cercana(C_Matrix::ElementT f);
void crear_matriz_nomax_orientacion(int i, int j);
void crear_matriz_nomax(int i, int j);
void seguir_cadena_m2(int i, int j);
void seguir_cadena_orientacion(int i, int j);
void juntar_contornos(int i, int j);

int tam, metodo, mascara, k;	//Tama�o, m�todo a utilizar, m�scara a utilizar, valor auxiliar para rellenar la m�scara 
double sigma;	//Valor de sigma (desviaci�n est�ndar)
int u_max, u_min;	//Umbrales m�ximos y m�nimos
C_Image imagen;	//Imagen original
C_Matrix kernel_gauss;	//Kernel de Gauss
C_Matrix matriz_J;	//Matriz convolucionada
C_Matrix matriz_es;	//Matriz de magnitud de los bordes
C_Matrix matriz_eo;	//Matriz de orientaci�n de los bordes
C_Matrix mascara_filtro_x;	//M�scara gen�rica X
C_Matrix mascara_filtro_y;	//M�scara gen�rica Y
C_Matrix matriz_Jx;	//Matriz convolucionada con m�scara X
C_Matrix matriz_Jy;	//M�triz convolucionada con m�scara Y
C_Matrix matriz_direccion;	//Matriz de direcci�n de los bordes
C_Matrix matriz_nomax;	//Matriz de no m�ximos
C_Matrix matriz_visitados;	//Matriz de p�xeles visitados
C_Matrix matriz_umbral;	//Matriz binaria resultante de la umbralizaci�n
//unsigned t0, t1;	//Variables medida de tiempo
C_Image imagen_final;	//Imagen final

int main(int argc, char **argv)
{
	printf("Inicio del programa\n");
	char txt_entrada[100];
	char txt_salida[100];
	cout << "Nombre de la imagen de entrada:\n";
	cin.getline(txt_entrada, 100);
	cout << "Nombre de la imagen de salida:\n";
	cin.getline(txt_salida, 100);
	if (strcmp(txt_entrada, txt_salida) == 0)printf("La imagen ser\240 sobrescrita\n");
	printf("Valor de sigma: ");
	scanf_s("%lf", &sigma);
	printf("Tama\244o kernel convoluci\242n (valor impar menor que la imagen): ");
	scanf_s("%d", &tam);
	while (true) {
		printf("Valor de umbral m\240ximo: ");
		scanf_s("%d", &u_max);
		printf("Valor de umbral m\241nimo: ");
		scanf_s("%d", &u_min);
		if (u_max > u_min) break;
		else printf("El valor m\240ximo debe ser mayor\n");
	}

	while (true) {
		cout << "Seleccionar m\240scara: 1 = Prewitt ; 2 = Sobel \n";
		scanf_s("%d", &mascara);
		if (mascara == 1 || mascara == 2) break;
		else printf("N\243mero 1 or 2\n");
	}
	while (true) {
		cout << "Seleccionar rutina: 1 = Orientaci\242n ; 2 = comprobar 8 vecinos:\n";
		scanf_s("%d", &metodo);
		if (metodo == 1 || metodo == 2) break;
		else printf("N\243mero 1 or 2\n");
	}

	//Medir el tiempo
	//t0 = clock();

	//Inicializaci�n del programa
	programa(txt_entrada, txt_salida);

	printf("Finalizado con \202xito\n");

	//Finalizaci�n medici�n del tiempo
	//t1 = clock();

	//double time = (double(t1 - t0) / CLOCKS_PER_SEC);
	//cout << "Tiempo de ejecuci\242n: " << time << endl;

}

/*	Programa principal que realiza el procesamiento de la imagen
	entrada[] : nombre de la imagen a procesar
	salida[] : nombre de la imagen creada
*/
void programa(char entrada[], char salida[]) {

	//Leer imagen de entrada
	printf("Leyendo imagen\n");

	imagen.ReadBMP(entrada);
	if (imagen.Fail()) {
		cout << "La imagen no existe.\n";
		cin.getline(entrada, 100);
		system("pause");
		exit(-1);
	}

	//Inicializar matrices
	kernel_gauss.Resize(1, tam, 1, tam);
	matriz_J.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	matriz_es.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	matriz_eo.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	mascara_filtro_x.Resize(1, 3, 1, 3);
	mascara_filtro_y.Resize(1, 3, 1, 3);
	matriz_Jx.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	matriz_Jy.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	matriz_direccion.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol());
	matriz_nomax.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);
	matriz_visitados.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol());
	matriz_umbral.Resize(imagen.FirstRow(), imagen.LastRow(), imagen.FirstCol(), imagen.LastCol(), 255);

	//Calcular kernel gaussiano
	//G(x,y) = (1/(2*PI*(sigma^2)))*e^-((x^2+y^2)/2+sigma^2)
	//sigma = desviaci�n est�ndar

	for (int i = kernel_gauss.FirstRow(); i <= kernel_gauss.LastRow(); i++) {
		for (int j = kernel_gauss.FirstCol(); j <= kernel_gauss.LastCol(); j++) {
			kernel_gauss(i, j) = (1 / ((2 * M_PI)*pow(sigma, 2))) * exp(-((pow(i, 2) + pow(j, 2)) / (2 * pow(sigma, 2))));
		}
	}

	kernel_gauss.DivideEscalar(kernel_gauss.Sum());

	//Matriz_J = I * G siendo I la original y G la gaussiana creada anteriormente
	//Hacer convoluci�n de I usando la m�scara G y guardarlo en J
	matriz_J = convolucion(imagen, kernel_gauss);

	//Gradiente Jx y Jy
	//Gradiente x e y de la imagen suavizada con el kernel gaussiano.

	/*
	M�scara 1 = Prewitt
	M�scara 2 = Sobel
	Default = Sobel
	*/

	k = mascara;
	if (k != 1 || k != 2) k = 2;

	//GradienteX
	mascara_filtro_x(1, 1) = -1;		mascara_filtro_x(1, 2) = 0;		mascara_filtro_x(1, 3) = 1;
	mascara_filtro_x(2, 1) = -k;		mascara_filtro_x(2, 2) = 0;		mascara_filtro_x(2, 3) = k;
	mascara_filtro_x(3, 1) = -1;		mascara_filtro_x(3, 2) = 0;		mascara_filtro_x(3, 3) = 1;

	//GradienteY
	mascara_filtro_y(1, 1) = -1;		mascara_filtro_y(1, 2) = -k;	mascara_filtro_y(1, 3) = -1;
	mascara_filtro_y(2, 1) = 0;			mascara_filtro_y(2, 2) = 0;		mascara_filtro_y(2, 3) = 0;
	mascara_filtro_y(3, 1) = 1;			mascara_filtro_y(3, 2) = k;		mascara_filtro_y(3, 3) = 1;

	//Convoluci�n Jx
	//Hacer convoluci�n de J usando la m�scara seg�n el filtro elegido y guardarlo en matriz_Jx

	printf("Convoluci\242n con el gradiente X\n");
	matriz_Jx = convolucion(matriz_J, mascara_filtro_x);

	//Convoluci�n Jy
	//Hacer convoluci�n de J usando la m�scara seg�n el filtro elegido y guardarlo en matriz_Jy
	printf("Convoluci\242n con el gradiente Y\n");
	matriz_Jy = convolucion(matriz_J, mascara_filtro_y);

	printf("Calculando magnitud de los bordes\n");
	//Calcular la magnitud de los bordes 
	for (int i = matriz_es.FirstRow(); i <= matriz_es.LastRow(); i++) {
		for (int j = matriz_es.FirstCol(); j <= matriz_es.LastCol(); j++) {
			matriz_es(i, j) = sqrt((pow(matriz_Jx(i, j), 2) + (pow(matriz_Jy(i, j), 2))));
		}
	}

	///M�todo1
	if (metodo == 1) {

		printf("Estimando orientaci\242n de los bordes\n");
		//Estimar la orientaci�n de la normal de los bordes
		for (int i = matriz_eo.FirstRow(); i <= matriz_eo.LastRow(); i++) {
			for (int j = matriz_eo.FirstCol(); j <= matriz_eo.LastCol(); j++) {
				matriz_eo(i, j) = atan(matriz_Jy(i, j) / matriz_Jx(i, j));
			}
		}

		//Matriz de direcci�n
		//Estimar direcci�n posible entre 0 - 45 - 90 - 135
		printf("Estimando \240ngulo de los bordes\n");
		for (int i = matriz_eo.FirstRow(); i <= matriz_eo.LastRow(); i++) {
			for (int j = matriz_eo.FirstCol(); j <= matriz_eo.LastCol(); j++) {
				matriz_direccion(i, j) = direccion_cercana(matriz_eo(i, j));
			}
		}

		//Matriz no_max vecino m�s cercano
		printf("Calculando matriz de no m\240ximos\n");
		matriz_nomax.SetValue(0);

		for (int i = matriz_es.FirstRow(); i <= matriz_es.LastRow(); i++) {
			for (int j = matriz_es.FirstCol(); j <= matriz_es.LastCol(); j++) {
				crear_matriz_nomax_orientacion(i, j);
			}
		}

		//Hist�resis seg�n el umbral
		printf("Recorrer imagen a\244adiendo valores seg\243n el umbral\n");
		matriz_visitados.SetValue(0);
		matriz_umbral.SetValue(0);
		for (int i = imagen.FirstRow(); i <= imagen.LastRow(); i++) {
			for (int j = imagen.FirstCol(); j <= imagen.LastCol(); j++) {
				//Si se ha visitado el punto continua la ejecuci�n
				if (matriz_visitados(i, j) == 1)continue;
				//Los bordes se dejan igual para evitar problemas al visitar vecinos que est�n fuera de rango
				if (i != imagen.FirstRow() || i != imagen.LastRow() || j != imagen.LastCol() || j != imagen.FirstCol()) {
					if (matriz_nomax(i, j) >= u_max) seguir_cadena_orientacion(i, j);
				}
			}
		}

		//Uni�n de bordes
		for (int i = matriz_nomax.FirstRow(); i <= matriz_nomax.LastRow(); i++) {
			for (int j = matriz_nomax.FirstCol(); j <= matriz_nomax.LastCol(); j++) {
				if (matriz_nomax(i, j) >= u_max)juntar_contornos(i, j);
			}
		}
	}
	///M�todo2
	else if (metodo == 2) {

		//Comprobar el peso con los vecino en la direcci�n correcta
		printf("Calculando matriz de no m\240ximos\n");
		matriz_nomax.SetValue(0);

		for (int i = matriz_es.FirstRow(); i <= matriz_es.LastRow(); i++) {
			for (int j = matriz_es.FirstCol(); j <= matriz_es.LastCol(); j++) {
				crear_matriz_nomax(i, j);
			}
		}

		//Umbral hist�resis
		//Hist�resis
		printf("Hist\202resis una vez definido el umbral\n");
		matriz_visitados.SetValue(0);
		matriz_umbral.SetValue(0);
		for (int i = imagen.FirstRow(); i <= imagen.LastRow(); i++) {
			for (int j = imagen.FirstCol(); j <= imagen.LastCol(); j++) {
				//Si se ha visitado el punto continua la ejecuci�n
				if (matriz_visitados(i, j) == 1)continue;
				//Los bordes se dejan igual para evitar problemas al visitar vecinos que est�n fuera de rango
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

/*	M�todo para realizar convoluci�n sobre una matriz
	m1 : matriz sobre la que se realiza la convoluci�n
	m2 : m�scara kernel
	return aux : Matriz tras la convoluci�n
*/
C_Matrix convolucion(C_Matrix m1, C_Matrix m2) {

	int center = m2.RowN() / 2;	//N�mero de p�xeles de margen
	C_Matrix aux(m1.FirstRow(), m1.LastRow(), m1.FirstCol(), m1.LastCol(), 255);
	C_Matrix::ElementT sumatoria_convolucion;

	//Recorrer matriz imagen original
	for (int i = m1.FirstRow(); i <= m1.LastRow(); i++) {
		for (int j = m1.FirstCol(); j <= m1.LastCol(); j++) {
			sumatoria_convolucion = 0;
			//Si el punto se encuentra dentro del margen se deja el p�xel con el valor actual y se continua la ejecuci�n
			if (i <= center || i > (m1.LastRow() - center) || j <= center || j > (m1.LastCol() - center)) {
				sumatoria_convolucion = m1(i, j);
			}
			//En caso contrario se realiza la convoluci�n
			else {
				for (int k = m2.FirstRow(); k <= m2.LastRow(); k++) {
					for (int l = m2.FirstCol(); l <= m2.LastCol(); l++) {
						//Se resta center por el tama�o del kernel
						sumatoria_convolucion += m1(i + k - center - m2.FirstRow(), j + l - center - m2.FirstCol())*m2(k, l);
					}
				}
			}
			//Se a�ade el valor calculado al punto correspondiente
			aux(i, j) = sumatoria_convolucion;
		}
	}

	return aux;
}

/*	M�todo para comprobar a qu� �ngulo se acerca m�s de 0�, 45�, 90� o 135�
	f : valor a estudiar
	return : �ngulo al que m�s se aproxima la entrada
*/
int direccion_cercana(C_Matrix::ElementT f) {

	//Convertir valor en �ngulo
	C_Matrix::ElementT angulo = (f / M_PI) * 180.0;
	//Comprobar cercan�a
	if ((angulo < 22.5 && angulo > -22.5) || (angulo > 157.5 && angulo < -157.5)) return 0;
	if ((angulo > 22.5 && angulo < 67.5) || (angulo < -112.5 && angulo > -157.5)) return 45;
	if ((angulo > 67.5 && angulo < 112.5) || (angulo < -67.5 && angulo > -112.5)) return 90;
	if ((angulo > 112.5 && angulo < 157.5) || (angulo < -22.5 && angulo > -67.5)) return 135;

	return -1;	//No llega aqu�
}

/*	M�todo para suprimir m�ximos por interpolaci�n de pesos seg�n cuadrante.
	i : posici�n i (filas) de la matriz
	j : posici�n j (columnas) de la matriz
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
			//Cuadrante 2 y 4
			if (matriz_Jx(i, j)*matriz_Jy(i, j) > 0) {
				//Izquierda arriba
				v1 = matriz_es(i - 1, j - 1);
				//Derecha abajo
				v3 = matriz_es(i + 1, j + 1);
			}
			//Cuadrante 1 y 3
			else {
				//Derecha arriba
				v1 = matriz_es(i - 1, j + 1);
				//Izquierda abajo
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
			//Izquierda y derecha
			v2 = matriz_es(i, j + 1);
			v4 = matriz_es(i, j - 1);
			//Cuadrante 3 y 1
			if (matriz_Jx(i, j)*matriz_Jy(i, j) > 0) {
				//Derecha abajo
				v1 = matriz_es(i + 1, j + 1);
				//Izquierda arriba
				v3 = matriz_es(i - 1, j - 1);
			}
			//Cuadrante 2 y 4
			else {
				//Derecha arriba
				v1 = matriz_es(i - 1, j + 1);
				//Izquierda abajo
				v3 = matriz_es(i + 1, j - 1);
			}
		}
		//Comparaci�n de p�xeles
		aux1 = peso*v1 + (1 - peso)*v2;
		aux2 = peso*v3 + (1 - peso)*v4;

		//Comparaci�n del m�ximo con el local
		//Igualar
		if (aux >= aux1 && aux >= aux2) {
			matriz_nomax(i, j) = aux;
		}
		//Suprimir
		else if (aux < aux1 || aux < aux2) {
			matriz_nomax(i, j) = 0;
		}
	}
}

/*	M�todo para crear la imagen seg�n el umbral de m�nimos de forma recursiva siguiendo los 8 vecinos
	i : posici�n fila de la matriz
	j : posici�n columna de la matriz
*/
void seguir_cadena_m2(int i, int j) {

	C_Matrix::ElementT valor;

	if (matriz_visitados(i, j) == 0)matriz_visitados(i, j) = 1;
	else return;

	//Marcar como borde
	matriz_umbral(i, j) = 255;

	//Vecinos del punto X
	int pos_x[8] = { 1,1,0,-1,-1,-1,0,1 };
	int pos_y[8] = { 0,1,1,1,0,-1,-1,-1 };

	//Se recorren los vecinos sin importar la direcci�n del punto
	for (int k = 0; k < 8; k++) {

		valor = matriz_nomax(i + pos_x[k], j + pos_y[k]);
		//Los que superen el umbral m�nimo pertenecen al borde tambi�n
		if (valor >= u_min) {
			seguir_cadena_m2(i + pos_x[k], j + pos_y[k]);
		}
	}
}

/*	M�todo para suprimir m�ximos dependiendo de la direcci�n
	i : posici�n i (filas) de la matriz
	j : posici�n j (columnas) de la matriz
*/
void crear_matriz_nomax_orientacion(int i, int j) {

	//Evitar p�xeles cercanos al borde para no salirse del rango de la imagen al bucar sus vecinos
	if (i == matriz_es.FirstRow() || i == matriz_es.LastRow() || j == matriz_es.FirstCol() || j == matriz_es.LastCol()) {
		matriz_nomax(i, j) = 0;
		return;
	}

	int direccion = matriz_direccion(i, j);
	switch (direccion) {
	case 0:	//Comprobar con los p�xeles de la izquierda y la derecha
		if (matriz_es(i, j) < matriz_es(i, j - 1) || matriz_es(i, j) < matriz_es(i, j + 1)) {
			matriz_nomax(i, j) = 0;
		}
		else {
			matriz_nomax(i, j) = matriz_es(i, j);
		}
		break;
	case 45://Comprobar con los p�xeles de la izquierda abajo y la derecha arriba
		if (matriz_es(i, j) < matriz_es(i - 1, j + 1) || matriz_es(i, j) < matriz_es(i + 1, j - 1)) {
			matriz_nomax(i, j) = 0;
		}
		else {
			matriz_nomax(i, j) = matriz_es(i, j);
		}
		break;
	case 90://Comprobar con los p�xeles de arriba y abajo
		if (matriz_es(i, j) < matriz_es(i - 1, j) || matriz_es(i, j) < matriz_es(i + 1, j)) {
			matriz_nomax(i, j) = 0;
		}
		else {
			matriz_nomax(i, j) = matriz_es(i, j);
		}
		break;
	case 135://Comprobar con los p�xeles de la izquierda arriba y la derecha abajo
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

/*	M�todo para crear la imagen seg�n el umbral de m�nimos de forma recursiva seg�n la orientaci�n
	i : posici�n fila de la matriz
	j : posici�n columna de la matriz
*/
void seguir_cadena_orientacion(int i, int j) {

	if (matriz_visitados(i, j) == 1)return;	//P�xel ya estudiado
	if (i == imagen.FirstRow() || i == imagen.LastRow() || j == imagen.FirstCol() || j == imagen.LastCol()) return;

	matriz_visitados(i, j) == 1;	//Visitado
	matriz_umbral(i, j) = 255;		//Marcado como borde

	//A partir de aqu� recorrer p�xeles conectados en ambas direcciones perpendiculares
	//a la normal del borde mientras sea > u_min

	//Valores de los dos vecinos
	int aux_x1, aux_y1, aux_x2, aux_y2;

	int direccion = matriz_direccion(i, j);
	switch (direccion) {
	case 0:	//Comprobar con los p�xeles de la izquierda y la derecha
		aux_x1 = 0; aux_x2 = 0; aux_y1 = -1; aux_y2 = 1;
		break;
	case 45://Comprobar con los p�xeles de la izquierda abajo y la derecha arriba
		aux_x1 = -1; aux_x2 = 1; aux_y1 = -1; aux_y2 = -1;
		break;
	case 90://Comprobar con los p�xeles de arriba y abajo
		aux_x1 = -1; aux_x2 = 1; aux_y1 = 0; aux_y2 = 0;
		break;
	case 135://Comprobar con los p�xeles de la izquierda arriba y la derecha abajo
		aux_x1 = -1; aux_x2 = 1; aux_y1 = -1; aux_y2 = 1;
		break;
	default:
		printf("Error\n");
		aux_x1 = 0; aux_x2 = 0; aux_y1 = 0; aux_y2 = 0;
		break;
	}
	//Seguir cadena por los puntos donde el valor sea mayor al umbral m�nimo
	if (matriz_nomax(i + aux_x1, j + aux_y1) >= u_min) seguir_cadena_orientacion(i + aux_x1, j + aux_y1);
	if (matriz_nomax(i + aux_x2, j + aux_y2) >= u_min) seguir_cadena_orientacion(i + aux_x2, j + aux_y2);
}

/*	M�todo para unir contornos d�biles pr�ximos a un contorno fuerte
	i: posici�n fila de la matriz
	j: posici�n columna de la matriz
*/
void juntar_contornos(int i, int j) {

	if (i == matriz_nomax.FirstRow() || i == matriz_nomax.LastRow()
		|| j == matriz_nomax.FirstCol() || j == matriz_nomax.LastCol()) return;

	//Recorrer imagen original con una m�scara 3x3 y comprobar si hay alg�n borde fuerte
	for (int k = -1; k <= 1; k++) {
		for (int l = -1; l <= 1; l++) {
			if (k == 0 && l == 0) continue;
			if (matriz_nomax(i + k, j + l) >= u_min) matriz_umbral(i, j) = 255;
		}
	}
}