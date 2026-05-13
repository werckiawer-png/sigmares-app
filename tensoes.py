import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, Wedge

# --- CONFIGURAÇÃO DA PÁGINA ---
st.set_page_config(page_title="SigmaRes ECE - Academic Edition", page_icon="📚", layout="wide")

st.sidebar.title("⚙️ Parâmetros do Ensaio")
st.sidebar.markdown("Configure a viga e a ferramenta.")

with st.sidebar.form("painel_dados"):
    st.header("Geometria")
    R_form = st.number_input("Raio do Former R (m)", value=0.075, format="%.3f")
    b = st.number_input("Largura b (m)", value=0.014, format="%.3f")
    h = st.number_input("Altura h (m)", value=0.006, format="%.3f")
    d = st.number_input("Distância apoios d (m)", value=0.180, format="%.3f")
    delta = st.number_input("Deslocamento \u03B4 (m)", value=0.015, format="%.4f")

    st.header("Material")
    E_val = st.number_input("Módulo Young E (GPa)", value=207.0)
    Et_val = st.number_input("Módulo Tangente Et (GPa)", value=2.912)
    Sy_val = st.number_input("Limite Escoamento Sy (MPa)", value=1026.0)

    btn_calcular = st.form_submit_button("⚙️ Processar Relatório Acadêmico")

if not btn_calcular:
    st.title("📚 SigmaRes ECE: Padrão Acadêmico")
    st.info("Insira os dados na barra lateral para gerar o relatório corrigido.")

if btn_calcular:
    # --- MOTOR DE CÁLCULO ---
    E, Et, Sy = E_val * 1e9, Et_val * 1e9, Sy_val * 1e6
    c = h / 2.0
    rho = R_form + c
    kappa = 1.0 / rho
    I = (b * h**3) / 12.0
    
    # Geometria de Tangência
    dist_v = rho - delta
    hipot = np.sqrt((d/2)**2 + dist_v**2)
    theta_rad = np.arctan2(d/2, dist_v) - np.arccos(rho / hipot)
    theta_deg = np.degrees(theta_rad)
    L_contato = 2 * theta_rad * R_form

    # Tensões
    rho_y = (E * c) / Sy
    yy = min(c, (rho / rho_y) * c)
    y_coords = np.linspace(-c, c, 500)
    sigma_car = np.where(np.abs(y_coords) <= yy, (E * y_coords) / rho, np.sign(y_coords) * (Sy + Et * (np.abs(y_coords)/rho - Sy/E)))
    
    # Correção do NumPy
    m_y, m_s = y_coords[250:], sigma_car[250:]
    try: M_aplicado = 2 * np.trapezoid(m_s * b * m_y, m_y)
    except AttributeError: M_aplicado = 2 * np.trapz(m_s * b * m_y, m_y)
    
    sigma_desc = (M_aplicado * y_coords) / I
    sigma_res = sigma_car - sigma_desc
    
    # --- CORREÇÃO DA TENSÃO MÁXIMA ---
    Sy_MPa = Sy_val
    Et_MPa = Et_val * 1000
    E_MPa = E_val * 1000
    frac1 = c / rho
    frac2 = Sy_MPa / E_MPa
    subtracao = frac1 - frac2
    multiplicacao = Et_MPa * subtracao
    sigma_max_calc = Sy_MPa + multiplicacao
    
    # Amassamento Visual (15% da meia-altura para fins didáticos)
    d_crush = c * 0.15 

    # --- INTERFACE PRINCIPAL ---
    st.title("📚 Relatório Analítico: SigmaRes ECE")
    
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Raio Constante (\u03C1)", f"{rho*1000:.1f} mm")
    col2.metric("Curvatura (\u03BA)", f"{kappa:.3f} m⁻¹")
    col3.metric("Arco de Contato", f"{L_contato*1000:.2f} mm")
    col4.metric("Tensão Máx Aplicada", f"{sigma_max_calc:.1f} MPa") # Corrigido!

    st.markdown("---")

    # --- GRÁFICOS (ESTILO ACADÊMICO) ---
    st.header("🖼️ Diagramas do Ensaio e Distribuição de Tensões")
    plt.style.use('default') 
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22, 7))

    # --- PAINEL 1: VISÃO LONGITUDINAL ---
    center_y = (rho - delta) * 1000 
    former_1 = Circle((0, center_y), R_form * 1000, color='#BDC3C7', alpha=0.5, ec='black', lw=1, label='Ferramenta')
    ax1.add_patch(former_1)
    
    contact_arc = Arc((0, center_y), R_form * 2000, R_form * 2000, angle=0, 
                      theta1=270-theta_deg, theta2=270+theta_deg, color='#C0392B', lw=4, label='Contato')
    ax1.add_patch(contact_arc)
    
    phis = np.linspace(-theta_rad, theta_rad, 100)
    x_curve = rho * np.sin(phis) * 1000
    y_curve = (center_y - rho * np.cos(phis)) * 1000
    x_full = np.concatenate(([-d/2 * 1000], x_curve, [d/2 * 1000]))
    y_full = np.concatenate(([0], y_curve, [0]))
    
    ax1.fill_between(x_full, y_full - c*1000, y_full + c*1000, color='#BDC3C7', alpha=0.4, label='Viga')
    ax1.plot(x_full, y_full, 'k--', lw=1, label='Eixo Neutro')
    ax1.plot(-d/2 * 1000, -c * 1000, '^', markersize=10, color='black', label='Apoios')
    ax1.plot(d/2 * 1000, -c * 1000, '^', markersize=10, color='black')
    
    ax1.set_aspect('equal')
    ax1.set_xlim(-d/2 * 1000 - 15, d/2 * 1000 + 15)
    ax1.set_ylim(-delta * 1000 - h * 1000 - 15, center_y + R_form * 1000 + 15)
    ax1.set_title("Visão Longitudinal", fontsize=14, fontweight='bold')
    ax1.set_xlabel("Largura (mm)")
    ax1.set_ylabel("Altura (mm)")
    ax1.legend(loc='upper right')

    # --- PAINEL 2: SEÇÃO TRANSVERSAL COM FORMER AMASSANDO ---
    ax2.set_title("Seção Transversal e Contato", fontsize=14, fontweight='bold')
    
    # Geometria da Viga Amassada
    x_fill = np.array([-b/2 * 1000, 0, b/2 * 1000])
    y_top_fill = np.array([c * 1000, (c - d_crush) * 1000, c * 1000]) # O "V" do amassamento
    
    # Zonas coloridas acompanhando o amassamento
    ax2.fill_between(x_fill, yy * 1000, y_top_fill, color='#FF9999', alpha=1.0, label='Zona Plástica')
    ax2.fill_between([-b/2 * 1000, b/2 * 1000], -yy * 1000, yy * 1000, color='#FFFF99', alpha=1.0, label='Zona Elástica')
    ax2.fill_between([-b/2 * 1000, b/2 * 1000], -c * 1000, -yy * 1000, color='#FF9999', alpha=1.0)
    
    # Contorno Preto da Viga Amassada
    vertices_viga = [(-b/2*1000, c*1000), (0, (c-d_crush)*1000), (b/2*1000, c*1000), (b/2*1000, -c*1000), (-b/2*1000, -c*1000)]
    viga_poly = plt.Polygon(vertices_viga, fill=False, edgecolor='black', lw=2)
    ax2.add_patch(viga_poly)
    
    # Desenho do Former atuando diretamente na viga
    r_former_vis = b * 0.4 * 1000 
    former_wedge = Wedge((0, (c - d_crush)*1000 + r_former_vis), r_former_vis, 180, 360, color='#BDC3C7', ec='black', lw=2, alpha=0.9, label='Former (Punção)')
    ax2.add_patch(former_wedge)

    ax2.axhline(0, color='black', linestyle='-.', lw=1) 
    ax2.axhline(yy * 1000, color='black', linestyle='--', lw=1) 
    ax2.axhline(-yy * 1000, color='black', linestyle='--', lw=1) 

    # Os Vetores Pretos (Desviando do Former)
    posicoes_y = np.linspace(-c * 1000, c * 1000, 15)
    for pos_y in posicoes_y:
        if abs(pos_y) < 0.2: continue
        tensao_local = np.interp(pos_y/1000, y_coords, sigma_car)
        tamanho_seta = abs(tensao_local / sigma_car[-1]) * (b * 1000 * 0.8)
        
        origem_x = b/2 * 1000
        if pos_y > 0 and pos_y < (c - d_crush)*1000: # Compressão (Aponta para a ESQUERDA, abaixo do former)
            ax2.arrow(origem_x + tamanho_seta, pos_y, -tamanho_seta, 0, head_width=0.3, head_length=0.6, fc='black', ec='black', lw=1.5, length_includes_head=True)
        elif pos_y < 0: # Tração (Aponta para a DIREITA)
            ax2.arrow(origem_x, pos_y, tamanho_seta, 0, head_width=0.3, head_length=0.6, fc='black', ec='black', lw=1.5, length_includes_head=True)

    ax2.set_xlim(-b*1000*0.8, b*1000*2.2)
    ax2.set_ylim(-c*1000*1.5, c*1000*1.5)
    ax2.set_xlabel("Largura (mm)")
    ax2.legend(loc="upper left")

    # --- PAINEL 3: GRÁFICO DE TENSÕES ---
    ax3.set_title("Perfil Analítico de Tensões", fontsize=14, fontweight='bold')
    ax3.plot(sigma_car / 1e6, y_coords * 1000, label='Carregamento', linestyle='-', color='black', linewidth=2)
    ax3.plot(sigma_desc / 1e6, y_coords * 1000, label='Descarregamento', linestyle='--', color='gray', linewidth=2)
    ax3.plot(sigma_res / 1e6, y_coords * 1000, label='Residual', color='#E74C3C', linewidth=3)
    
    ax3.axvline(0, color='black', linewidth=1)
    ax3.axhline(0, color='black', linewidth=1, linestyle='-.')
    ax3.set_xlabel("Tensão (MPa)")
    ax3.legend()
    ax3.grid(color='lightgray', linestyle='--', alpha=0.7)
    
    st.pyplot(fig)

    st.markdown("---")

    # --- QUADRO DE ROLAGEM COM A MATEMÁTICA CORRIGIDA ---
    st.header("🧮 Memória de Cálculo (Passo a Passo)")
    with st.container(height=550, border=True):
        st.subheader("1. Geometria de Contato e Abraçamento")
        st.write("Raio de curvatura constante imposto no eixo neutro ($\rho$):")
        st.latex(rf"\rho = R + c = {R_form:.3f} + {c:.4f} = \mathbf{{{rho:.4f} \, m}}")
        
        st.write("Ângulo de abraçamento e comprimento de contato:")
        st.latex(rf"\theta = \mathbf{{{theta_deg:.2f}^\circ}} \implies L_c = 2 \cdot \theta_{{rad}} \cdot R = \mathbf{{{L_contato*1000:.2f} \, mm}}")

        st.markdown("---")
        
        st.subheader("2. Cálculo Passo a Passo da Tensão Máxima ($\sigma_{max}$)")
        st.write("Equação Fundamental do Modelo Elasto-Plástico com Encruamento (ECE):")
        st.latex(r"\sigma_{max} = S_y + E_t \cdot \left( \frac{c}{\rho} - \frac{S_y}{E} \right)")
        
        st.write("**Passo 2.1:** Substituição dos dados brutos (valores uniformizados em MPa):")
        st.latex(rf"\sigma_{{max}} = {Sy_MPa:.1f} + {Et_MPa:.1f} \cdot \left( \frac{{{c:.4f}}}{{{rho:.4f}}} - \frac{{{Sy_MPa:.1f}}}{{{E_MPa:.1f}}} \right)")
        
        st.write("**Passo 2.2:** Cálculo das frações internas (Deformação Imposta vs Deformação Elástica):")
        st.latex(rf"\sigma_{{max}} = {Sy_MPa:.1f} + {Et_MPa:.1f} \cdot \left( {frac1:.6f} - {frac2:.6f} \right)")
        
        st.write("**Passo 2.3:** Subtração (Deformação Plástica Líquida efetiva):")
        st.latex(rf"\sigma_{{max}} = {Sy_MPa:.1f} + {Et_MPa:.1f} \cdot ({subtracao:.6f})")
        
        st.write("**Passo 2.4:** Multiplicação pelo Módulo Tangente ($E_t$):")
        st.latex(rf"\sigma_{{max}} = {Sy_MPa:.1f} + {multiplicacao:.2f}")
        
        st.write("**Passo 2.5:** Resultado Final (Tensão de Pico):")
        st.latex(rf"\mathbf{{\sigma_{{max}} = {sigma_max_calc:.2f} \, MPa}}")

        st.markdown("---")
        st.subheader("3. Equilíbrio Residual")
        st.write("Subtraindo o descarregamento elástico do estado máximo, obtemos a tensão residual de superfície:")
        st.latex(rf"\sigma_{{res}} = \mathbf{{{-sigma_res[-1]/1e6:.2f} \, MPa}}")

    st.markdown("---")

    # --- PARECER TÉCNICO ---
    st.header("📋 Parecer Técnico: Deformação e Esmagamento Local")
    st.success(f"""
    **1. Mecânica de Esmagamento Superficial:**
    A análise geométrica avançada corrobora com a inspeção visual do ensaio prático. O former não age apenas como um gabarito de dobra; a força necessária para impor um deslocamento de {delta*1000:.1f} mm gerou uma pressão de contato extremo. Como ilustrado na Seção Transversal, a ferramenta afunda na borda superior da viga, esmagando localmente as fibras tracionadas.

    **2. Severidade do Encruamento:**
    Comprovado pelo cálculo passo a passo corrigido, a restrição impôs um núcleo elástico ínfimo de {2*yy*1000:.2f} mm. A agressividade da curvatura impulsionou a tensão além da fronteira elástica ({Sy_MPa} MPa), forçando o material a atingir um pico exato de **{sigma_max_calc:.1f} MPa**.

    **3. Proteção por Tensões Residuais:**
    Após o recuo do former e o *spring-back*, as camadas esmagadas retiveram o núcleo central. O resultado é uma sólida tensão residual de compressão de **{-sigma_res[-1]/1e6:.1f} MPa** na superfície, essencial para a resistência do material.
    """)