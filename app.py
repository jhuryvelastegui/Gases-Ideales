import streamlit as st
import matplotlib.pyplot as plt
import time
import requests

# Configuración de la página web
st.set_page_config(
    page_title="Leyes de los Gases Ideales - Estudio",
    layout="centered",
    initial_sidebar_state="expanded"
)

R = 0.082  # L·atm/(mol·K)

# Inicializar variables de estado (Session State de Streamlit)
if "tiempo_inicio" not in st.session_state:
    st.session_state.tiempo_inicio = None
if "ejercicio_en_curso" not in st.session_state:
    st.session_state.ejercicio_en_curso = False
if "resultado_texto" not in st.session_state:
    st.session_state.resultado_texto = ""
if "grafica_args" not in st.session_state:
    st.session_state.grafica_args = None

def linspace(a, b, n):
    if n < 2: return [a]
    return [a + (b - a) * i / (n - 1) for i in range(n)]

# ================== DEFINICIÓN DE LEYES ==================
LEYES = {
    "Boyle (P₁V₁ = P₂V₂)": {
        "vars": ["p1","v1","p2","v2"],
        "labels": {"p1":"P₁ (atm)","v1":"V₁ (L)","p2":"P₂ (atm)","v2":"V₂ (L)"},
        "calculos": {
            "p1": lambda d: (d["p2"]*d["v2"])/d["v1"],
            "v1": lambda d: (d["p2"]*d["v2"])/d["p1"],
            "p2": lambda d: (d["p1"]*d["v1"])/d["v2"],
            "v2": lambda d: (d["p1"]*d["v1"])/d["p2"],
        },
        "unidades": {"p1":"atm","v1":"L","p2":"atm","v2":"L"},
        "ecuacion": "P₁V₁ = P₂V₂",
        "grafica": lambda d, obj, res: graficar_boyle(d, obj, res)
    },
    "Charles (V₁/T₁ = V₂/T₂)": {
        "vars": ["v1","t1","v2","t2"],
        "labels": {"v1":"V₁ (L)","t1":"T₁ (K)","v2":"V₂ (L)","t2":"T₂ (K)"},
        "calculos": {
            "v1": lambda d: (d["v2"]*d["t1"])/d["t2"],
            "t1": lambda d: (d["v1"]*d["t2"])/d["v2"],
            "v2": lambda d: (d["v1"]*d["t2"])/d["t1"],
            "t2": lambda d: (d["v2"]*d["t1"])/d["v1"],
        },
        "unidades": {"v1":"L","t1":"K","v2":"L","t2":"K"},
        "ecuacion": "V₁/T₁ = V₂/T₂",
        "grafica": lambda d, obj, res: graficar_charles(d, obj, res)
    },
    "Gay-Lussac (P₁/T₁ = P₂/T₂)": {
        "vars": ["p1","t1","p2","t2"],
        "labels": {"p1":"P₁ (atm)","t1":"T₁ (K)","p2":"P₂ (atm)","t2":"T₂ (K)"},
        "calculos": {
            "p1": lambda d: (d["p2"]*d["t1"])/d["t2"],
            "t1": lambda d: (d["p1"]*d["t2"])/d["p2"],
            "p2": lambda d: (d["p1"]*d["t2"])/d["t1"],
            "t2": lambda d: (d["p2"]*d["t1"])/d["p1"],
        },
        "unidades": {"p1":"atm","t1":"K","p2":"atm","t2":"K"},
        "ecuacion": "P₁/T₁ = P₂/T₂",
        "grafica": lambda d, obj, res: graficar_gaylussac(d, obj, res)
    },
    "Combinada (P₁V₁/T₁ = P₂V₂/T₂)": {
        "vars": ["p1","v1","t1","p2","v2","t2"],
        "labels": {"p1":"P₁ (atm)","v1":"V₁ (L)","t1":"T₁ (K)","p2":"P₂ (atm)","v2":"V₂ (L)","t2":"T₂ (K)"},
        "calculos": {
            "p1": lambda d: (d["p2"]*d["v2"]*d["t1"])/(d["v1"]*d["t2"]),
            "v1": lambda d: (d["p2"]*d["v2"]*d["t1"])/(d["p1"]*d["t2"]),
            "t1": lambda d: (d["p1"]*d["v1"]*d["t2"])/(d["p2"]*d["v2"]),
            "p2": lambda d: (d["p1"]*d["v1"]*d["t2"])/(d["v2"]*d["t1"]),
            "v2": lambda d: (d["p1"]*d["v1"]*d["t2"])/(d["p2"]*d["t1"]),
            "t2": lambda d: (d["p2"]*d["v2"]*d["t1"])/(d["p1"]*d["v1"]),
        },
        "unidades": {"p1":"atm","v1":"L","t1":"K","p2":"atm","v2":"L","t2":"K"},
        "ecuacion": "P₁V₁/T₁ = P₂V₂/T₂",
        "grafica": lambda d, obj, res: graficar_combinada(d, obj, res)
    },
    "Gas ideal (PV = nRT)": {
        "vars": ["p","v","n","t"],
        "labels": {"p":"P (atm)","v":"V (L)","n":"n (mol)","t":"T (K)"},
        "calculos": {
            "p": lambda d: (d["n"]*R*d["t"])/d["v"],
            "v": lambda d: (d["n"]*R*d["t"])/d["p"],
            "n": lambda d: (d["p"]*d["v"])/(R*d["t"]),
            "t": lambda d: (d["p"]*d["v"])/(d["n"]*R),
        },
        "unidades": {"p":"atm","v":"L","n":"mol","t":"K"},
        "ecuacion": "PV = nRT",
        "grafica": lambda d, obj, res: graficar_gasideal(d, obj, res)
    },
}

# ================== FUNCIONES DE GRAFICACIÓN ==================
def graficar_boyle(d, objetivo, resultado):
    p1,v1,p2,v2 = d["p1"],d["v1"],d["p2"],d["v2"]
    k = p1*v1
    v_min = max(0.05, min(v1,v2)*0.3)
    v_max = max(v1,v2)*2.2
    Vs = linspace(v_min, v_max, 120)
    return [
        {"x": Vs, "y": [k/v for v in Vs], "color": "#1e40af", "label": "Isoterma"},
        {"x": [v1], "y": [p1], "color": "#dc2626", "label": "Estado 1", "points_only": True},
        {"x": [v2], "y": [p2], "color": "#16a34a", "label": "Estado 2", "points_only": True},
    ], "Volumen (L)", "Presión (atm)", "Boyle: P ∝ 1/V  (T cte)"

def graficar_charles(d, objetivo, resultado):
    v1,t1,v2,t2 = d["v1"],d["t1"],d["v2"],d["t2"]
    k = v1/t1
    t_min = max(1, min(t1,t2)*0.5)
    t_max = max(t1,t2)*1.6
    Ts = linspace(t_min, t_max, 120)
    return [
        {"x": Ts, "y": [k*t for t in Ts], "color": "#1e40af", "label": "Isobara"},
        {"x": [t1], "y": [v1], "color": "#dc2626", "label": "Estado 1", "points_only": True},
        {"x": [t2], "y": [v2], "color": "#16a34a", "label": "Estado 2", "points_only": True},
    ], "Temperatura (K)", "Volumen (L)", "Charles: V ∝ T  (P cte)"

def graficar_gaylussac(d, objetivo, resultado):
    p1,t1,p2,t2 = d["p1"],d["t1"],d["p2"],d["t2"]
    k = p1/t1
    t_min = max(1, min(t1,t2)*0.5)
    t_max = max(t1,t2)*1.6
    Ts = linspace(t_min, t_max, 120)
    return [
        {"x": Ts, "y": [k*t for t in Ts], "color": "#1e40af", "label": "Isocora"},
        {"x": [t1], "y": [p1], "color": "#dc2626", "label": "Estado 1", "points_only": True},
        {"x": [t2], "y": [p2], "color": "#16a34a", "label": "Estado 2", "points_only": True},
    ], "Temperatura (K)", "Presión (atm)", "Gay-Lussac: P ∝ T  (V cte)"

def graficar_combinada(d, objetivo, resultado):
    p1,v1,t1 = d["p1"],d["v1"],d["t1"]
    p2,v2,t2 = d["p2"],d["v2"],d["t2"]
    v_min = max(0.05, min(v1,v2)*0.3)
    v_max = max(v1,v2)*2.5
    Vs = linspace(v_min, v_max, 120)
    k1 = p1*v1; k2 = p2*v2
    return [
        {"x": Vs, "y": [k1/v for v in Vs], "color": "#dc2626", "label": f"Isoterma T₁={t1:.0f}K"},
        {"x": Vs, "y": [k2/v for v in Vs], "color": "#16a34a", "label": f"Isoterma T₂={t2:.1f}K"},
        {"x": [v1], "y": [p1], "color": "#dc2626", "label": "Estado 1", "points_only": True},
        {"x": [v2], "y": [p2], "color": "#16a34a", "label": "Estado 2", "points_only": True},
    ], "Volumen (L)", "Presión (atm)", "Ley combinada — diagrama P-V"

def graficar_gasideal(d, objetivo, resultado):
    p,v,n,t = d["p"],d["v"],d["n"],d["t"]
    if objetivo in ("p","v"):
        v_min = max(0.01, v*0.2)
        v_max = v*3
        Vs = linspace(v_min, v_max, 120)
        return [
            {"x": Vs, "y": [(n*R*t)/vv for vv in Vs], "color": "#1e40af", "label": f"T={t:.1f}K, n={n:.3g}mol"},
            {"x": [v], "y": [p], "color": "#16a34a", "label": "Punto actual", "points_only": True},
        ], "Volumen (L)", "Presión (atm)", "Gas ideal: P ∝ 1/V"
    elif objetivo == "t":
        p_min = max(0.01, p*0.2)
        p_max = p*3
        Ps = linspace(p_min, p_max, 120)
        return [
            {"x": Ps, "y": [(pp*v)/(n*R) for pp in Ps], "color": "#1e40af", "label": f"V={v:.3g}L, n={n:.3g}mol"},
            {"x": [p], "y": [t], "color": "#16a34a", "label": "Punto actual", "points_only": True},
        ], "Presión (atm)", "Temperatura (K)", "Gas ideal: T ∝ P"
    else:
        t_min = max(1, t*0.3)
        t_max = t*2.5
        Ts = linspace(t_min, t_max, 120)
        return [
            {"x": Ts, "y": [(p*v)/(R*tt) for tt in Ts], "color": "#1e40af", "label": f"P={p:.3g}atm, V={v:.3g}L"},
            {"x": [t], "y": [n], "color": "#16a34a", "label": "Punto actual", "points_only": True},
        ], "Temperatura (K)", "Cantidad (mol)", "Gas ideal: n ∝ 1/T"

# Función puente para renderizar con Matplotlib en la Web
def renderizar_grafica_web(data, xl, yl, title):
    fig, ax = plt.subplots(figsize=(7, 3.8))
    for s in data:
        if s.get("points_only"):
            ax.scatter(s["x"], s["y"], color=s["color"], label=s["label"], s=100, edgecolors="white", zorder=5)
        else:
            ax.plot(s["x"], s["y"], color=s["color"], label=s["label"], linewidth=2)
    ax.set_title(title, fontsize=11, fontweight="bold", color="#111827")
    ax.set_xlabel(xl, fontsize=9, color="#6b7280")
    ax.set_ylabel(yl, fontsize=9, color="#6b7280")
    ax.grid(True, linestyle="--", alpha=0.5, color="#d1d5db")
    ax.legend(fontsize=8)
    
    st.pyplot(fig)
    plt.close(fig)  # Liberar memoria

# ================== INTERFAZ DE USUARIO (STREAMLIT) ==================

st.title("🧪 Leyes de los Gases Ideales")
st.subheader("Módulo de Estudio e Interacción Laboratorial")

# --- Barra Lateral ---
with st.sidebar:
    st.markdown("### Instrucciones")
    st.info("Haga click en el botón rojo, Iniciar Ejercicio, para desbloquear la celdas de ingreso de datos y comenzar la resolución de su ejercicio.")
    
    st.divider()
    st.markdown("### ⏱️ Control de Cronómetro")
    
    # Inyección de CSS para asegurar que el botón de Iniciar sea rojo
    st.markdown("""
    <style>
    div[data-testid="stButton"] button[kind="primary"] {
        background-color: #dc2626 !important;
        border-color: #dc2626 !important;
        color: white !important;
    }
    </style>
    """, unsafe_allow_html=True)
    
    # Lógica de botones del cronómetro
    if not st.session_state.ejercicio_en_curso:
        if st.button("▶ Iniciar Ejercicio", type="primary", use_container_width=True):
            st.session_state.tiempo_inicio = time.time()
            st.session_state.ejercicio_en_curso = True
            st.session_state.resultado_texto = ""
            st.session_state.grafica_args = None
            st.rerun()
    else:
        st.success("⏱️ Ejercicio en curso...")
        if st.button("⏹️ Cancelar Tiempo", use_container_width=True):
            st.session_state.tiempo_inicio = None
            st.session_state.ejercicio_en_curso = False
            st.rerun()

    # --- Créditos Chiquitos (Footer Lateral) ---
    st.markdown("<br><br><br><br>", unsafe_allow_html=True) # Espacio para empujarlo abajo
    st.markdown("""
    <div style='font-size: 11px; color: #9ca3af; text-align: left; line-height: 1.4;'>
        J. Velasteguí<br>
        M. Quijia<br>
        D. Pila<br>
        Universidad Central del Ecuador
    </div>
    """, unsafe_allow_html=True)

# --- Configuración del problema ---
col_ley, col_var = st.columns(2)

with col_ley:
    ley_seleccionada = st.selectbox("Seleccione la Ley:", list(LEYES.keys()))
    ley = LEYES[ley_seleccionada]

with col_var:
    opciones_var = [f"{ley['labels'][v]} ← calcular" for v in ley["vars"]]
    var_idx = st.selectbox("Variable a calcular:", range(len(opciones_var)), format_func=lambda x: opciones_var[x])
    objetivo = ley["vars"][var_idx]

# Mostrar Ecuación Destacada
st.info(f"**Ecuación base:** ` {ley['ecuacion']} `", icon="💡")

# --- Formulario de Datos de Entrada ---
st.markdown("#### Datos de entrada")
datos_ingresados = {}

columnas_formulario = st.columns(2)
for index, var in enumerate(ley["vars"]):
    col_destino = columnas_formulario[index % 2]
    with col_destino:
        if var == objetivo:
            st.text_input(f"{ley['labels'][var]} (Resultado)", value="Se calculará automáticamente...", disabled=True)
        else:
            deshabilitado = not st.session_state.ejercicio_en_curso
            val = st.number_input(f"Ingrese {ley['labels'][var]}:", value=0.0, step=0.1, min_value=0.0, disabled=deshabilitado, key=f"input_{ley_seleccionada}_{var}")
            datos_ingresados[var] = val if val != 0.0 else None

# --- Botón de Ejecución de Cálculo ---
if st.session_state.ejercicio_en_curso:
    if st.button("🧮 Calcular Resultado", use_container_width=True):
        faltan_datos = False
        for v in ley["vars"]:
            if v != objetivo and datos_ingresados.get(v) is None:
                st.error(f"❌ Falta ingresar el valor para: {ley['labels'][v]}")
                faltan_datos = True
        
        if not faltan_datos:
            try:
                calc_dict = {k: v for k, v in datos_ingresados.items() if v is not None}
                resultado = ley["calculos"][objetivo](calc_dict)
                
                # Calcular Tiempos
                tiempo_finalizacion = time.time()
                tiempo_total = tiempo_finalizacion - st.session_state.tiempo_inicio
                respuesta_str = f"{resultado:.5g}"
                
                # ==========================================
                # Envío silencioso a Google Forms (6 variables)
                # ==========================================
                form_url = "https://docs.google.com/forms/d/e/1FAIpQLSelDUqWAKuPihWFDerdoq5_VwzYfwuJpXLrZIuSfK-WHDCpYA/formResponse"
                
                payload = {
                    "entry.913368875": "Anónimo",                 # A: Nombre (Ya no se pide)
                    "entry.2011290175": "N/A",                    # B: Carrera (Ya no se pide)
                    "entry.1519319718": ley_seleccionada,         # C: Ley Evaluada
                    "entry.1124505437": objetivo.upper(),         # D: Variable Calculada
                    "entry.812643603": f"{tiempo_total:.2f}",     # E: Tiempo en segundos
                    "entry.790733520": respuesta_str              # F: Respuesta Obtenida
                }
                
                try:
                    requests.post(form_url, data=payload, timeout=5)
                except Exception:
                    pass
                # ==========================================
                
                # Guardar estados para mostrar en pantalla
                unidad = ley["unidades"][objetivo]
                st.session_state.resultado_texto = f"✅ **Resultado:** {ley['labels'][objetivo].split('(')[0].strip()} = {respuesta_str} {unidad}  |  ⏱️ **Tiempo total:** {tiempo_total:.2f} s"
                
                datos_completos = {**calc_dict, objetivo: resultado}
                st.session_state.grafica_args = ley["grafica"](datos_completos, objetivo, resultado)
                
                st.session_state.ejercicio_en_curso = False
                st.session_state.tiempo_inicio = None
                st.rerun()
                
            except ZeroDivisionError:
                st.error("❌ Error matemático: División por cero.")
            except Exception as e:
                st.error(f"❌ Ocurrió un error en el cálculo: {e}")
else:
    if st.session_state.resultado_texto == "":
        st.caption("ℹ️ Por favor, presione 'Iniciar Ejercicio' en la barra lateral izquierda para desbloquear los campos.")

# --- Sección de Resultados y Gráficas ---
if st.session_state.resultado_texto:
    st.success(st.session_state.resultado_texto)

if st.session_state.grafica_args:
    with st.expander("📊 Ver Representación Gráfica", expanded=True):
        data, xl, yl, title = st.session_state.grafica_args
        renderizar_grafica_web(data, xl, yl, title)
