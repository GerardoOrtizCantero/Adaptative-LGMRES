# PD-LGMRES: Variante Adaptativa del LGMRES

Este repositorio contiene una implementación en MATLAB del método **PD-LGMRES**, una extensión adaptativa del LGMRES(m,k) que ajusta automáticamente el parámetro de reinicio `m` usando un controlador Proporcional‑Derivativo (PD).

## ¿Qué contiene?

- `pd_lgmres.m`: Implementación del método PD‑LGMRES
- `pd_rule.m`: Función que aplica la regla PD para actualizar `m`
- Ejemplos de uso y scripts para replicar los resultados del trabajo de tesis

## Antecedentes

El PD‑LGMRES se creo usando como referencia los solvers alojados en el **KrySBAS** (Krylov Subspace-Based Adaptive Solvers), mantenida por el **Núcleo de Investigación y Desarrollo Tecnológico (NIDTEC)** de la Facultad Politécnica, Universidad Nacional de Asunción, Paraguay. 

- Repositorio oficial: [nidtec-una/krysbas-dev](https://github.com/nidtec-una/krysbas-dev)  

## Instalación

```sh
git clone https://github.com/tu-usuario/pd-lgmres.git
cd pd-lgmres
addpath('pd-lgmres')
